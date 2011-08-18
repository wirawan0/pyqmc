# $Id: slaterdet.py,v 1.4 2011-08-18 16:17:36 wirawan Exp $
#
# pyqmc.matrices.slaterdet
# Created: 20110617
# Wirawan Purwanto
# College of William and Mary
# Williamsburg, VA, USA
#
# This module is part of PyQMC project.
#

"""
pyqmc.matrices.slaterdet

Module for handling a generic Slater determinant.
"""

import numpy
from numpy import abs as Abs
# scipy is too sweeping, we will revert to numpy instead:
from numpy import array, asmatrix, asarray, conj, matrix, ndarray, sqrt
from numpy import product
from numpy.linalg import det

class Det(object):
  """Basic Slater determinant object.
  Object members:
  - up = up determinant: matrix
  - dn = down determinant: matrix
  """
  def __init__(self, src_up=None, src_dn=None, ndn=None):
    if src_up != None:
      self.make(src_up, src_dn, ndn)
  def norm(self):
    return (det(self.up.H * self.up) \
            * det(self.dn.H * self.dn)) # .real
  def make(self, src_up, src_dn = None, ndn = None):
    self.up = asmatrix(array(src_up,order='C'))
    if (src_dn != None):
      self.dn = asmatrix(array(src_dn,order='C'))
    elif (ndn != None):
      self.dn = self.up[:, 0:ndn]
    else:
      self.dn = self.up
    if (self.up.shape[0] != self.dn.shape[0]):
      raise ValueError("Det: mismatch number of basis funcs for up and down dets")
  def nup(self):
    return self.up.shape[1]
  def ndn(self):
    return self.dn.shape[1]
  def nbasis(self):
    return self.up.shape[0]

  # provide alias for some newer APIs that require this:
  @property
  def alpha(self):
    return self.up
  @alpha.setter
  def alpha(self, x):
    self.up = x
  @property
  def beta(self):
    return self.dn
  @beta.setter
  def beta(self, x):
    self.dn = x


class MultiDet(object):
  """Basic multi-Slater determinant object.
  Object members:
  * dets[:] = a list of Det objects, with attached "ampl"
  # field in it.
  """
  def set_ampl(self,mtx):
    '''Sets the determinant amplitudes. The mtx argument is either a numpy
    array, or python list or tuple.'''
    if isinstance(mtx, ndarray):
      for (d, ampl) in zip(self.dets, asarray(mtx).reshape(len(self.dets))):
        d.ampl = ampl
    elif isinstance(mtx, (list, tuple)):
      for (d, ampl) in zip(self.dets, mtx):
        d.ampl = ampl
  def make_ci(self, ampl, cfg, orbs, nup, ndn):
    '''Creates a multideterminant wfn by using CI-like configuration.
    Arguments:
    - ampl: a list of amplitudes, either as 1-dimensional array/matrix or
      Python list.
    - cfg: list of determinant configurations. This should be a 2D array
      where each row is a configuration of the up and down electrons.
    - orbs: a list of single-particle orbitals, given as a 2D array
      (where each column is a single-particle orbital).
    '''
    ndets = len(cfg)
    # this will enforce ampl to have the same length as ndets (and also
    # make it insensitive to orientation (whether ndets-length array,
    # or (1 x ndets) or (ndets x 1) matrix.
    ampl = asarray(ampl).reshape(ndets)
    orbs = asarray(orbs)
    nptot = nup + ndn
    nptot2 = len(cfg[0])
    if nptot != nptot2:
      raise ValueError, \
        "Mismatch total number of particles in `cfg' to that specified in nup+ndn"

    up_cols = xrange(0, nup)
    dn_cols = xrange(nup, nup+ndn)
    up_shape = (orbs.shape[0], nup)
    dn_shape = (orbs.shape[0], ndn)

    self.dets = []
    #asarray([ orbs.T[i] for i in c[:nup] ]).T,
    #asarray([ orbs.T[j] for j in c[nup:] ]).T
    # I use for loop here so we know if there's an error in the
    # multidet configuration:
    for c in cfg:
      newdet = Det(
                   asarray([ orbs.T[i] for i in c[:nup] ]).T,
                   asarray([ orbs.T[j] for j in c[nup:] ]).T
                  )
      self.dets.append(newdet)

    self.set_ampl(ampl)
    self.cfg_ci = asarray(cfg)
    self.ampl_ci = asarray(ampl)
  def norm(self):
    return mdet_ovlp(self)
  def normalize(self, N = None):
    if (N != None):
      n = Abs(N)
    else:
      n = mdet_ovlp(self)
    self.norm_orig = n
    for d in self.dets: d.ampl *= 1.0 / sqrt(n)
  @property
  def ndets(self):
    return len(self.dets)

  # Allow access in array-like fashion:
  def __len__(self):
    return len(self.dets)
  def __getitem__(self, i):
    return self.dets[i]
  def __iter__(self):
    return self.dets.__iter__()

# Functions:

def up_ovlp(det1, det2):
  '''Computes overlap between two determinants (up dets only)'''
  return det( det1.up.H * det2.up )

def dn_ovlp(det1, det2):
  '''Computes overlap between two determinants (down dets only)'''
  return det( det1.dn.H * det2.dn )

def det_ovlp(det1, det2):
  '''Computes overlap between two determinants (Det objects)'''
  return det( det1.up.H * det2.up ) \
       * det( det1.dn.H * det2.dn )

def mdet_ovlp(mdet1, mdet2 = None):
  '''Computes overlap between two MultiDet objects'''
  ovlp = 0.0 + 0.0j
  if (mdet2 == None):
    L = 0
    for Ldet in mdet1.dets:
      ovlp += abs(Ldet.ampl)**2 * Ldet.norm()
      L += 1
      for Rdet in mdet1.dets[L:]:
        ovlp += 2 * (conj(Ldet.ampl) * Rdet.ampl \
                     * det_ovlp(Ldet, Rdet)).real
  else:
    for Ldet in mdet1.dets:
      for Rdet in mdet2.dets:
        ovlp += conj(Ldet.ampl) * Rdet.ampl \
              * det_ovlp(Ldet, Rdet)
  return ovlp

