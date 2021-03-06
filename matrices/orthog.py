# $Id: orthog.py,v 1.2 2011-06-20 16:17:52 wirawan Exp $
#
# pyqmc.matrices.orthog
# Created: 20110617
# Wirawan Purwanto
# College of William and Mary
# Williamsburg, VA, USA
#
# Basis/orbital orthogonalization module.
# This module is part of PyQMC project.
#

"""
pyqmc.matrices.orthog

Module for handling basis/orbital/vector orthogonalization.
"""

import numpy
from numpy import sqrt
from numpy import diag, dot
from numpy.linalg import eigh
from wpylib.math.linalg.gram_schmidt import modgs

matmul = numpy.dot

def Xorth_canonical(S):
  """Makes the canonical orthogonalizer.
  Returns a tuple containing the X matrix and its inverse.
  """
  (s, U) = eigh(S)
  sqrt_s1 = sqrt(s)
  inv_sqrt_s1 = 1.0 / sqrt_s1
  X = matmul(U, diag(inv_sqrt_s1))
  Xinv = matmul(diag(sqrt_s1), U.T.conj())
  return (X, Xinv)


def normalize_eigenvectors(U):
  """Normalizes a matrix of eigenvalues (stored in columns)
  so that the diagonal values are all positive."""
  from numpy import abs
  rslt = numpy.empty_like(U)
  for i in xrange(len(U)):
    Uii = U[i,i]
    if Uii != 0.0:
      U_norm = abs(Uii) / Uii
      rslt[:,i] = U[:,i] * U_norm
  return rslt

