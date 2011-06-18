# $Id: gms.py,v 1.1 2011-06-18 04:01:53 wirawan Exp $
#
# pyqmc.matrices.gms
# Created: 20110617
# Wirawan Purwanto
# College of William and Mary
# Williamsburg, VA, USA
#
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


