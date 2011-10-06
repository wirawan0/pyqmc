# $Id: __init__.py,v 1.1 2011-10-06 19:41:50 wirawan Exp $
#
# pyqmc.wfn.planewave module
# Created: 20100929
# Wirawan Purwanto
#

"""
Plane-wave type wave functions.
Mostly for reading wave functions produced by other computed packages:
ABINIT, PWSCF, etc.

"""

from wpylib.iofmt.text_input import text_input
from scipy.linalg.fblas import zgemm

class pw_text1(object):
  """Plane-wave orbitals.
  The file is a plain text file, containing rows of seven fields:

     counter   Gx  Gy  Gz  (nonblank,unused)  Real_ampl  Imag_ampl

  There are nbasis rows per orbital.
  """
  def __init__(self, nbasis):
    self.nbasis = nbasis

  def load_orbital(self, infile=None):
    """Loads an orbital from a planewave orbital list.
    The `infile' argument must be a text_input object."""
    if infile == None: infile = self.infile
    rslt = infile.read_items((1,int,'Gx'), (2,int,'Gy'), (3,int,'Gz'),
                             ((5,6),complex,'phi'), maxcount=self.nbasis)
    return rslt
