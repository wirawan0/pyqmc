# $Id: orthog.py,v 1.1 2011-06-18 04:00:37 wirawan Exp $
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
