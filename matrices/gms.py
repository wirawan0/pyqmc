# $Id: gms.py,v 1.14 2011-09-12 21:59:24 wirawan Exp $
#
# pyqmc.matrices.gms
# Created: 20110617
# Wirawan Purwanto
# College of William and Mary
# Williamsburg, VA, USA
#
# Wissam's "GAMESS"-style matrix handling.
# This module is part of PyQMC project.
#

"""
pyqmc.matrices.gms

Module for handling the so-called `GAMESS-style' matrices, i.e. the
legacy matrix format used by GAFQMC code.
"""

import sys
import os

import numpy
import numpy.linalg

from numpy import asmatrix, empty, matrix, zeros_like
from numpy.linalg import eigh, inv

import wpylib.math
from wpylib.sugar import ifelse
from wpylib.iofmt.fortbin import fortran_bin_file
from wpylib.iofmt.text_input import text_input
from wpylib.iofmt.text_output import text_output
from pyqmc.matrices.orthog import Xorth_canonical
from pyqmc.matrices.utils import read_matrix, read_det_matrix

from pyqmc import PyqmcDataError

class OneBodyGms(object):
  """One-body matrix element.
  This will read the `overlap' (S) and `core hamiltonian' (H1) matrices
  from a given input file.
  """

  def __init__(self, src=None):
    if src:
      self.read(src)

  def read(self, infile):
    from os.path import abspath
    F = text_input(infile)
    self.filename = infile
    self.filename_abs = abspath(infile)
    r = F.next_rec()
    self.nbasis = int(r[0])
    self.nelem = int(r[1])
    (n1, m1) = self.read_matrix(F)
    (n2, m2) = self.read_matrix(F)
    self.S = m1
    self.H1 = m2
    self.S_name = n1
    self.H1_name = n2
    F.close()

  def write(self, outfile, comment=None):
    F = text_output(outfile)
    if comment:
      cmt = " # " + str(comment)
    else:
      cmt = ""
    nbasis = self.nbasis
    nelem = nbasis * (nbasis + 1) / 2
    F.write("%d %d%s\n" % (nbasis, nelem, cmt))
    self.write_matrix(F, self.S, getattr(self, "S_name", "Overlap"), symmetric=True)
    self.write_matrix(F, self.H1, getattr(self, "H1_name", "Core Hamiltonian"), symmetric=True)
    F.close()

  def read_matrix(self, F):
    """Reads an individual matrix from text stream F, each line in the
    format:

        row col value
        row col value
        ...

    there are self.nelem such rows."""
    matname = F.next()
    m1 = F.read_items((0,int), (1,int), (2,float), maxcount=self.nelem)
    m2 = numpy.zeros((self.nbasis,self.nbasis), dtype=float)
    for (r,c,v) in m1:
      m2[r-1,c-1] = v
      m2[c-1,r-1] = v
    return (matname, m2)

  def write_matrix(self, F, m, name, symmetric=True):
    """Writes a matrix in one_body_gms format.
    If symmetric == True, then m[r,c] is written where r >= c is used.
    """
    m = numpy.asarray(m)
    r_len = len("%d" % m.shape[0])
    c_len = len("%d" % m.shape[1])
    fmt = "%" + str(r_len) + "d %" + str(c_len) + "d %-.15g\n"

    F.write("%s\n" % str(name))
    for (r,R) in enumerate(m):
      if symmetric:
        c_max = r+1
      else:
        c_max = len(R)
      for (c,RC) in enumerate(R[:c_max]):
        F.write(fmt % (r+1, c+1, RC))
    F.write("\n")
    F.flush()


  def compute_Xorth(self):
    (X,Xinv) = Xorth_canonical(self.S)
    self.Xorth_ = X
    self.Xorth_inv_ = Xinv

  @property
  def Xorth(self):
    if not hasattr(self, "Xorth_"):
      self.compute_Xorth()
    return self.Xorth_

  @property
  def Xorth_inv(self):
    if not hasattr(self, "Xorth_inv_"):
      self.compute_Xorth()
    return self.Xorth_inv_


class TwoBodyGmsUfmt(object):
  """Two-body matrix elements.
  WARNING: This piece of code is extremely slow and the matrix can be
  very large.
  The matrix member `H2' is a four-dimension object written in quantum
  Chemistry's notation:

  (il|jk)  ==> V_ijkl ==> <ij|lk>  ==> H2[i,l,j,k]

  """
  V2b_permutation_options = {
    0: 0,
    False: 0,
    None: 0,
    1: 1,
    'standard': 1,
    'general': 1,
    2: 2,
    'real': 2,
    'default': 2, # the default option!
  }

  def init_record_info(self):
    """Numpy datatype representing a two_body_gms_ufmt record entry, which
    *includes* the Fortran record marker."""
    lsb_msb = '<' # LSB default
    int32 = lsb_msb + 'i4'
    float64 = lsb_msb + 'f8'
    self.rectype = numpy.dtype([('m1', int32), \
        ('i', int32), ('l', int32), ('j', int32), ('k', int32), \
        ('v', float64), \
        ('m2', int32)])
    # net bytesize excluding record markers
    self.recsize_net = self.rectype.itemsize - numpy.dtype(int32).itemsize * 2
  def make_flat_index(self, i,l,j,k):
    """Make flat indices.
    For a 4-dim array containing self.nbasis elements in each dimension.
    WARNING: C ordering is assumed here."""
    nb = self.nbasis
    return ((i * nb + l) * nb + j) * nb + k
  def __init__(self, infile=None, nbasis=None, debug=None, blksize=16384, perm='default'):
    self.init_record_info()
    if infile:
      self.read(infile, nbasis, debug=debug, blksize=blksize, perm=perm)
  def read(self, infile, nbasis, debug=None, blksize=16384, perm='default'):
    """Reads in the matrix elements from a Fortran binary file.
    This is supposed to be an accelerated implementation.
    We *bypass* the fortran binary format and slurp the file into memory
    before doing further processing.

    Permutation flags (`perm`) honored:
    * 0, False, None = No permutation  (generally you don't want this except for
      debugging)
    * 1, 'standard' = standard fourfold permutation for a Hermitian two-body
      Hamiltonian
    * 2, 'real' = eightfold permutation for a Hermitian two-body
      Hamiltonian with real basis functions in real space
    """
    from numpy import conj
    from os.path import abspath
    assert nbasis > 0
    try:
      perm = self.V2b_permutation_options[perm]
    except KeyError:
      raise ValueError, "Invalid permutation options: `%s'" % perm
    self.nbasis = nbasis
    dbg = text_output(ifelse(debug, sys.stdout, None), flush=True)
    H2 = numpy.zeros((nbasis,nbasis,nbasis,nbasis), dtype=float)
    self.H2 = H2
    nn = nbasis * (nbasis + 1) // 2
    S = os.stat(infile)
    fsize = S.st_size
    # net bytesize excluding marker
    rec_count = fsize // self.rectype.itemsize
    dbg("File %s: %d integral records to be read\n" % (infile, rec_count))
    dbg("Matrix element permutation flag = %s\n" % (perm))
    F = open(infile, "rb")
    self.filename = infile
    self.filename_abs = abspath(infile)
    # We use blocked read and assignment to minimize the python overhead
    for iblk in xrange(0, rec_count, blksize):
      read_blksize = min(blksize, rec_count - iblk)
      blob = numpy.fromfile(F, dtype=self.rectype, count=read_blksize)
      # The following provides a minimal consistency check if the file just read
      # is indeed a valid two_body_gms_ufmt file:
      if not numpy.all(blob['m1'] == self.recsize_net):
        raise PyqmcDataError, \
          "Invalid record marker (m1) detected: file %s may be corrupt or of incorrect format." \
          % (infile,)
      if not numpy.all(blob['m2'] == self.recsize_net):
        raise PyqmcDataError, \
          "Invalid record marker (m2) detected: file %s may be corrupt or of incorrect format." \
          % (infile,)
      # convert to py index (0-based)
      blob['i'] -= 1
      blob['l'] -= 1
      blob['j'] -= 1
      blob['k'] -= 1
      get_flat_perm_index = lambda iljk: \
        self.make_flat_index(blob[iljk[0]], blob[iljk[1]], blob[iljk[2]], blob[iljk[3]])
      # Use: V2b_inspect.permute_V2b('i','j','l','k',chem=1)
      # to generate this:
      v = blob['v']
      H2.put(get_flat_perm_index('iljk'), v)
      if (perm == 1 or perm == 2):
        H2.put(get_flat_perm_index('jkil'), v)
        H2.put(get_flat_perm_index('likj'), conj(v))
        H2.put(get_flat_perm_index('kjli'), conj(v))
        if (perm == 2):
          # Only usable if the basis orbitals are real in real-space.
          H2.put(get_flat_perm_index('ilkj'), v)
          H2.put(get_flat_perm_index('kjil'), v)
          H2.put(get_flat_perm_index('jkli'), v)
          H2.put(get_flat_perm_index('lijk'), v)
    F.close()
  def read0(self, infile, nbasis, debug=None):
    """Reads in the matrix elements from a Fortran binary file.
    This is a reference implementation, very slow.
    """
    assert nbasis > 0
    dbg = text_output(ifelse(debug, sys.stdout, None), flush=True)
    H2 = numpy.zeros((nbasis,nbasis,nbasis,nbasis), dtype=float)
    self.H2 = H2
    F = fortran_bin_file(infile)
    S = os.stat(infile)
    fsize = S.st_size
    rec_desc = (('i', numpy.int32, 4), ('v', numpy.float64))
    rec_bytesize = F.byte_length(*rec_desc) + F.default_int(0).itemsize * 2
    rec_count = fsize // rec_bytesize
    dbg("File %s: %d integral records to be read\n" % (infile, rec_count))
    for cc in xrange(rec_count):
      rec = F.read(*rec_desc)
      rec['i'] -= 1 # convert to py index
      (i,l,j,k) = rec['i']
      v = rec['v']
      # Use: V2b_inspect.permute_V2b('i','j','l','k',chem=1)
      # to generate this:
      H2[i,l, j,k] = v
      H2[j,k, i,l] = v
      H2[l,i, k,j] = v
      H2[k,j, l,i] = v
      H2[i,l, k,j] = v
      H2[k,j, i,l] = v
      H2[j,k, l,i] = v
      H2[l,i, j,k] = v
    F.close()
  def write(self, outfile, eps=wpylib.math.epsilon(numpy.float64)):
    """Writes the matrix elements in the two_body_gms_ufmt binary format.
    The eigthfold symmetry of the integral is considered here, and values
    smaller than double precision `eps' are ignored.
    """
    from numpy import abs
    H2 = self.H2
    nbasis = H2.shape[0]
    F = fortran_bin_file(outfile, "w")
    # Must preserve the dtype of i,l,j,k; so we have to use arange
    # and not do arithmetic on i,j,k,l when writing them out:
    ibasis = numpy.arange(1, nbasis+1, dtype=numpy.int32)
    # In LD format: (i,l) is taken only if (i => l)
    for k in ibasis:
      for j in ibasis[k-1:]:
        LD_jk = LD(j-1,k-1,nbasis)
        # LD(j,k)
        #if j < k: continue
        for l in ibasis:
          for i in ibasis[l-1:]:
            #if i < l: continue
            LD_il = LD(i-1,l-1,nbasis)
            if LD_il < LD_jk: continue
            v = H2[i-1,l-1, j-1,k-1]
            if abs(v) < eps: continue
            F.write_vals(i,l, j,k, v)
    F.close()


class EigenGms(object): #{
  """Molecular orbitals expressed in (nonorthogonal) AO basis---the so-called
  `eigen_gms' format.
  """
  def __init__(self, src=None):
    if src:
      self.read(src)

  def read(self, infile, fort70=None, verbose=None): #{
    '''Reads in eigen_gms formatted file.
    The parsing results are stored in:
    rslt.alpha A (nbasis x norb) matrix for alpha orbitals in (nonorthogonal)
               AO basis
    rslt.beta  A (nbasis x norb) matrix for beta orbitals in (nonorthogonal)
               AO basis
    rslt.udet  If nonzero, we have just read separate alpha and beta orbitals.
               Otherwise, beta field is not defined.
    BE WARNED: the parser is very primitive.
    It's not intended to be foolproof or flexible.
    '''
    from os.path import abspath
    if hasattr(infile, "next"):
      self_open = False
    elif isinstance(infile, basestring): # if a string, let's open the file
      if verbose: print "Reading eigen_gms from file " + infile
      self.filename = infile
      self.filename_abs = abspath(infile)
      infile = open(infile, "r")
      self_open = True
    else:
      raise TypeError, \
        "The infile argument must be either an open file, " + \
        "a stream with next() method, or a string (filename)"

    # The first text line must be correct:
    fld = infile.next().split()
    self.nbasis = int(fld[0])
    self.norb = int(fld[1])
    if verbose:
      print "EigenGms.read: nbasis=%d, norb=%d" % (self.nbasis, self.norb)

    infile.next() # skip a blank line
    self.alpha = asmatrix(empty((self.nbasis, self.norb)))
    for o in xrange(0, self.norb):
      self.alpha[:,o] = \
        read_matrix(infile, "eigen_gms_alpha_orb#" + str(o),
                    self.nbasis, 1)
      infile.next() # skip a blank line

    self.udet = False
    try:
      fld = infile.next().split()
      nbasis2 = int(fld[0])
      norb2 = int(fld[1])
      if (nbasis2 != self.nbasis or norb2 != self.norb):
        sys.stderr.write("Warning: invalid number of basis funcs " + \
                         "or orbitals for beta sector; " + \
                         "omitting beta sector altogether.")
      else:
        infile.next() # skip a blank line
        self.beta = asmatrix(empty((self.nbasis, self.norb)))
        for o in xrange(0, self.norb):
          self.beta[:,o] = \
            read_matrix(infile, "eigen_gms_beta_orb#" + str(o),
                        self.nbasis, 1)
          infile.next() # skip a blank line
        self.udet = True
    except:
      # If it doesn't find the beta orbital, forget about it.
      pass

    if self_open: infile.close()

    if verbose:
      if (self.udet):
        print "Both alpha and beta orbital sectors were read"
      else:
        print "Only alpha orbital sector was read"

    if fort70 != None:
      if verbose:
        print "fort70 provided, will do orthogonalization"
      iXorth = inv(fort70.Xorth())
      self.alpha = iXorth * self.alpha;
      if hasattr(self, "beta"):
        self.beta = iXorth * self.beta;


  def write(self, outfile, comment=None, udet=False, verbose=None): #{
    """Writes orbitals in eigen_gms formatted file.

    If only `alpha' orbitals exist then only the alpha sector is written out.
    This can be overriden by setting udet==True; then the alpha sector is
    duplicated as beta as well.
    The `udet' argument is not used if `beta' exists; both sectors will
    always be written out.
    """
    out = text_output(outfile)

    (nbasis, norb) = self.alpha.shape
    if comment:
      cmt = " # " + str(comment)
    else:
      cmt = ""

    if hasattr(self, "beta"):
      sectors = (self.alpha, self.beta)
      if verbose:
        print "EigenGms.write: (alpha,beta) nbasis=%d, norb=%d" % (nbasis, norb)
    elif udet:
      sectors = (self.alpha, self.alpha)
      if verbose:
        print "EigenGms.write: (alpha,alpha) nbasis=%d, norb=%d" % (nbasis, norb)
    else:
      sectors = (self.alpha,)
      if verbose:
        print "EigenGms.write: (alpha only) nbasis=%d, norb=%d" % (nbasis, norb)

    for SS in sectors:
      out("%d %d%s\n\n" % (nbasis, norb, cmt))

      for orb in numpy.array(SS,copy=False).T:
        out("\n".join([ "%.15g" % a for a in orb ] + ["\n"]))



class Fort70(object): #{
  """Fort70 is a representation of the GAFQMC fort.70-style matrix file.
  As of 2011 this object is largely deprecated in favor of direct access
  to the source matrices via OneBodyGms and EigenGms objects.
  """
  #nup = None
  #ndn = None
  # Mapping between fort.70 matrix name and internal field name:
  matrix_map = {
    "S_eigen_val_sqrt" : "inv_sqrt_s",
    "S_eigen_vec" : "U",
    "H1" : "H1",
    "psiT_ampl" : "ampl",
    "MO_noise1_up": "noise_up",
    "MO_noise1_dn": "noise_dn",
  }
  matrix_complex = [
    "psiT_ampl",
  ]

  def __init__(self, fname = None, nup = None, ndn = None, verbose = None):
    if fname != None:
      self.nup = nup
      self.ndn = ndn
      self.read(fname, verbose)

  def check_nparts(self):
    if (self.nup == None or self.ndn == None):
      raise ValueError, \
        "Number of particles (nup, ndn) has not been initialized properly."
    #if (self.ndn > self.nup):
    #  raise ValueError, "The case of (ndn > nup) is not supported."

  def read(self, fname, verbose = None): #{
    '''Reads in fort.70-formatted file as returned by GAFQMC code.
    '''
    global psit
    #print "opening file " + fname
    self.check_nparts()
    nup = self.nup
    ndn = self.ndn

    # Always start afresh: delete all previously set data
    for k in Fort70.matrix_map.keys() + ["psiT_det", "Xorth_"]:
      if (hasattr(self, k)): delattr(self, k)

    #print "opening file " + fname
    inp = open(fname, "r")
    try:
      txt = inp.next().rstrip()
      if (txt != "GAFQMC matrix element file v1"):
        raise PyqmcDataError("Not a GAFQMC matrix element file: " + fname)
      txt = inp.next()
      # BE WARNED: the parser below is very primitive.
      # It's not intended to be foolproof.
      while (txt != ""):
        txt = txt.strip()
        if (txt.startswith("MATRIX ")):
          # new matrix is encountered
          fld = txt.split()
          name = fld[1]
          rows = int(fld[2])
          cols = int(fld[3])
          if verbose:
            print "Found matrix %s %d %d" % (name, rows, cols)
          if (name in Fort70.matrix_map):
            # Most matrices can be set using the mapping thing above:
            if name not in self.matrix_complex:
              setattr(self, self.matrix_map[name],
                      read_matrix(inp, name, rows, cols))
            else:
              mtx = read_matrix(inp, name, rows, cols*2)
              mtx = mtx[:,0::2] + 1j * mtx[:,1::2]
              setattr(self, self.matrix_map[name], mtx)
          elif (name.startswith("psiT_det_")):
            # The psiT number must be in order (1..NPsiTDet) or else the
            # code would read everything erroneously:
            det_no = int(name[9:]) - 1
            if (not hasattr(self, "psiT_det")): self.psiT_det = []
            if (cols == nup + ndn):
              restricted = False
            elif (cols == max(nup,ndn)):
              restricted = True
            else:
              raise PyqmcDataError, \
                ("Invalid number of columns for WF matrix `%s' (%dx%d): " + \
                "The valid ncols is either %d or %d") % \
                (name, rows, cols, nup+ndn,nup)
            detmp = read_det_matrix(inp, name, rows, nup, ndn,
                                    cplx=False, restricted=restricted,
                                    verbose=verbose)
            self.psiT_det.append(detmp)
          else:
            raise PyqmcDataError, "Unknown matrix name: `" + name + "'"
        txt = inp.next()
    except StopIteration:
      inp.close()
    if not inp.closed:
      inp.close()
  #}read

  @property
  def Xorth(self):
    if not hasattr(self, "Xorth_"):
      invs = matrix(zeros_like(self.U))
      for i in xrange(0, invs.shape[0]):
        invs[i,i] = self.inv_sqrt_s[i,0]
      self.Xorth_ = self.U * invs
      #print invs
      #for i in xrange(0, self.U.shape[1]):
      #  x = self.inv_sqrt_s[i,0]
      #  print i, x, type(x), type(self.Xorth_[:,i])
      #  self.Xorth_[:,i] *= self.inv_sqrt_s[i,0]
      #  print self.Xorth_[:,i]
    return self.Xorth_
#}Fort70 class


def LD(i,j,N):
  """python equivalent of gafqmc_LD on nwchem-gafqmc integral
  dumper module.
  Translates a lower-diagonal index (ii >= jj) to linear index
  0, 1, 2, 3, ...
  This follows python convention; thus 0 <= i < N, and so also j.
  This should be able to accept numpy.array datatype as i and j
  parameters.
  """
  # iskip is row traversal, jskip is column traversal.
  # (iskip+jskip) is the final array index.
  ii = numpy.maximum(i,j)
  jj = numpy.minimum(i,j)
  #if i >= j:
  #  ii = i
  #  jj = j
  #else:
  #  ii = j
  #  jj = i

  iskip = ii - jj # + 1
  #jskip = (jj-1)*N - (jj-2)*(jj-1)/2  # for 1-based
  jskip = (jj)*N - (jj-1)*(jj)//2  # for 0-based
  return iskip + jskip

