# $Id: utils.py,v 1.3 2011-06-18 04:12:14 wirawan Exp $
#
# pyqmc.matrices.utils
# Created: 20110617
# Wirawan Purwanto
# College of William and Mary
# Williamsburg, VA, USA
#
# This module is part of PyQMC project
#

"""
pyqmc.matrices.utils

General-purpose matrix-related utilities.
"""

import numpy
# I cross my finger: hopefully this is not circular:
from pyqmc.matrices.slaterdet import Det

## Matrix formation utilities

def complex_matrix(mat):
  """Converts real-imag interleaved real matrix into a python complex matrix."""
  return mat[:,0::2] + 1j * mat[:,1::2]


## Matrix reading utilities

def read_matrix(infile, name, rows, cols):
  """Reads a matrix from an opened text file.
  Currently, the `name' argument is used only to aid debugging (i.e. to
  give user more information concerning the matrix that fails to read).
  Changed: Now it uses `next()' method of the infile file-like object to
  extract text lines one by one.
  The matrix content must *not* contain any comment or stray field.
  The number of columns must be specified by the caller; but we can leave the
  number of rows unspecified, i.e. we read all the matrix contents until the
  end of the file.
  """
  r = 0
  rslt = []
  #print "read_matrix", name, rows, cols
  while rows == None or r < rows:
    try:
      txt = infile.next()
    except StopIteration:
      if (rows != None):
        raise EOFError, \
          "Unexpected EOF while reading the content of `%s' matrix (%sx%d): " + \
          "%d rows read so far" % (name, rows, cols, r)
      else:
        break
    # Remove comments -- UNDC: for later
    #txt = txt.split("!")[0]
    #txt = txt.split("#")[0]
    fld = txt.split()
    #print fld, len(fld), cols, (len(fld) - cols)
    if (len(fld) != cols):
      raise ValueError, \
        ("Incorrect number of columns for `%s' matrix (%sx%d): " + \
        "got %d columns at row %d") % \
        (name, rows, cols, len(fld), r+1)
    rslt.append(map(float, fld))
    r += 1
  rslt = numpy.matrix(rslt)
  return rslt


def read_det_matrix(infile, name, nbasis, nup, ndn, cplx=None, restricted=None, verbose=None):
  '''Reads a determinant from an opened text file, formatted as a matrix.
  This routine uses read_matrix to read in the matrix.
  '''
  cplx = bool(cplx)
  restricted = bool(restricted)
  if not restricted:
    cols = nup + ndn
  else:
    cols = max(nup,ndn)
  if cplx: cols *= 2
  if verbose:
    print "reading det matrix `%s' (%dx%d): %d basis, %d up, %d dn (%s, %s)" % \
          (name, nbasis, cols, nbasis, nup, ndn, \
           {True: "restricted", False: "unrestricted"}[restricted], \
           {True: "complex", False: "real"}[cplx])
  mtx = read_matrix(infile, name, nbasis, cols)
  if (cplx): mtx = mtx[:,0::2] + 1j * mtx[:,1::2]
  if (restricted):
    return Det(mtx[:, 0:nup], mtx[:, 0:ndn])
  else:
    return Det(mtx[:, 0:nup], mtx[:, nup:nup+ndn])


