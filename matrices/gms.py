# $Id: gms.py,v 1.2 2011-06-18 04:05:55 wirawan Exp $
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

import numpy
import numpy.linalg
from numpy.linalg import eigh

from wpylib.iofmt.text_input import text_input
from pyqmc.matrices.orthog import Xorth_canonical

class OneBodyGms(object):
  """One-body matrix element.
  This will get the `overlap' (S) and `core hamiltonian' (H1) matrices.
  """

  def __init__(self, src=None):
    if src:
      self.read(src)

  def read(self, infile):
    F = text_input(infile)
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


class EigenGms(object): #{
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
    if "readline" in dir(infile):
      self_open = False
    elif isinstance(infile, str): # if a string, let's open the file
      if verbose: print "Reading eigen_gms from file " + infile
      infile = open(infile, "r")
      self_open = True
    else:
      raise TypeError, \
        "The infile argument must be either an open file, " + \
        "a stream with readline() method, or a string (filename)"

    # The first text line must be correct:
    fld = infile.readline().split()
    self.nbasis = int(fld[0])
    self.norb = int(fld[1])
    if verbose:
      print "EigenGms.read: nbasis=%d, norb=%d" % (self.nbasis, self.norb)

    infile.readline() # skip a blank line
    self.alpha = asmatrix(empty((self.nbasis, self.norb)))
    for o in xrange(0, self.norb):
      self.alpha[:,o] = \
        read_matrix(infile, "eigen_gms_alpha_orb#" + str(o),
                    self.nbasis, 1)
      infile.readline() # skip a blank line

    self.udet = False
    try:
      fld = infile.readline().split()
      nbasis2 = int(fld[0])
      norb2 = int(fld[1])
      if (nbasis2 != self.nbasis or norb2 != self.norb):
        sys.stderr.write("Warning: invalid number of basis funcs " + \
                         "or orbitals for beta sector; " + \
                         "omitting beta sector altogether.")
      else:
        infile.readline() # skip a blank line
        self.beta = asmatrix(empty((self.nbasis, self.norb)))
        for o in xrange(0, self.norb):
          self.beta[:,o] = \
            read_matrix(infile, "eigen_gms_beta_orb#" + str(o),
                        self.nbasis, 1)
          infile.readline() # skip a blank line
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
  #} read_eigen_gms
#}EigenGms class


class Fort70(object): #{
  '''Fort70 is a representation of the GAFQMC fort.70-style matrix file.
  '''
  #nup = None
  #ndn = None
  # Mapping between fort.70 matrix name and internal field name:
  matrix_map = {
    "S_eigen_val_sqrt" : "inv_sqrt_s",
    "S_eigen_vec" : "U",
    "H1" : "H1",
    "psiT_ampl" : "ampl",
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
      txt = inp.readline().rstrip()
      if (txt != "GAFQMC matrix element file v1"):
        raise ValueError("Not a GAFQMC matrix element file: " + fname)
      txt = inp.readline()
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
              raise ValueError, \
                ("Invalid number of columns for WF matrix `%s' (%dx%d): " + \
                "The valid ncols is either %d or %d") % \
                (name, rows, cols, nup+ndn,nup)
            detmp = read_det_matrix(inp, name, rows, nup, ndn,
                                    cplx=False, restricted=restricted,
                                    verbose=verbose)
            self.psiT_det.append(detmp)
          else:
            raise ValueError, "Unknown matrix name: `" + name + "'"
        txt = inp.readline()
    finally:
      inp.close()
  #}read

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
