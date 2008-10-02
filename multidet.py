# $Id: multidet.py,v 1.1 2008-10-02 04:06:25 wirawan Exp $
#
# multidet.py
# Created: 20080925
# Wirawan Purwanto
#
# Determinant handling in Python.
# Unfortunately what looks nice in octave has its own cost: time.
# Hopefully this is not the case with python/scipy.

import scipy
import scipy.linalg
from numpy import matrix, conj
from numpy.linalg import det

class Det(object):
  # Contains:
  # up = up determinant: matrix
  # dn = down determinant: matrix
  def __init__(self, src_up = None, src_dn = None, ndn = None):
    if src_up != None:
      self.make(src_up, src_dn, ndn)
  def norm(self):
    return det( self.up.H * self.up ) \
         * det( self.dn.H * self.dn )
  def make(self, src_up, src_dn = None, ndn = None):
    self.up = matrix(src_up)
    if (src_dn != None):
      self.dn = matrix(src_dn)
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

class MultiDet(object):
  ampl = []
  det = []

# Functions:

def up_ovlp(det1, det2):
  '''Computes overlap between two determinants (up dets only)'''
  return det( det1.up.H * det2.up )

def dn_ovlp(det1, det2):
  '''Computes overlap between two determinants (dn dets only)'''
  return det( det1.dn.H * det2.dn )

def det_ovlp(det1, det2): #{
  '''Computes overlap between two determinants (Det objects)'''
  return det( det1.up.H * det2.up ) \
       * det( det1.dn.H * det2.dn )
#}det_ovlp

def mdet_ovlp(mdet1, mdet2 = None): #{
  '''Computes overlap between two MultiDet objects'''
  ovlp = 0 + 0j
  Lcnt = len(mdet1.ampl)
  if (mdet2 == None):
    L = 0
    for L in xrange(0, Lcnt):
      for R in xrange(L, Lcnt):
        ovlp1 = conj(mdet1.ampl[L]) * mdet1.ampl[R] \
              * det_ovlp(mdet1.det[L], mdet1.det[R])
        if (L == R):
          ovlp += ovlp1
        else:
          ovlp += conj(ovlp1) + ovlp1
  else:
    Rcnt = len(mdet2.ampl)
    for L in xrange(0, Lcnt):
      for R in xrange(0, Rcnt):
        ovlp += conj(mdet1.ampl[L]) * mdet1.ampl[R] \
              * det_ovlp(mdet1.det[L], mdet1.det[R])
  return ovlp
#}mdet_ovlp

def blah(det1):
  print det1

