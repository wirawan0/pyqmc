#
# pyqmc.matrices.dense module
# Created: 20140819
# Wirawan Purwanto
#
# Standard matrix handling.
# This module is part of PyQMC project.
#

"""
pyqmc.matrices.dense module
Dense matrix I/O handling.

"""



class mat2d(object):
  """A very simple Fortran binary format for a 2-D dense matrix.
  """
  debug = 0

  def read(self, F, **opts):
    from wpylib.iofmt.fortbin import fortran_bin_file, copy_subarray
    dbg = self.debug
    fbin = fortran_bin_file(F)
    self.fbin = fbin
    Int = fbin.default_int
    Dbl = float
    ndim = fbin.read(('ndim',Int))['ndim']
    if dbg >= 1: print "ndim = ", ndim
    if ndim != 2:
      raise PyqmcDataError, "Invalid array dimensionality."

    shape = fbin.read(('shape',Int,ndim))['shape']
    self.shape = shape
    if dbg >= 1: print "shape = ", shape
    maxcols = opts.get("max_cols", shape[1])
    if maxcols < shape[1]:
      shape = (shape[0], maxcols)
      if dbg >= 1: print "Limiting loaded data shape to ", shape
    dtype = self.guess_datatype_()
    # Define a fake datatype: an array of 'dtype' values, of length
    # 'shape[0]' for fast reading:
    dtype_M = numpy.dtype([('d', dtype, shape[0])])
    arr_tmp = fbin.bulk_read_array1(dtype=dtype_M, shape=(shape[1],))
    arr = copy_subarray(arr_tmp, 'd', order='F')
    return arr

  def guess_datatype_(self):
    """Internal subroutine to deduce what datatype is contained in
    the matrix file.
    """
    data_rec_len = self.fbin.peek_next_rec_len()
    if self.debug >= 1: print "data rec len = ", data_rec_len
    M = self.shape[0]
    if data_rec_len == numpy.float64(0).itemsize * M:
      dtype = numpy.float64
    elif data_rec_len == numpy.complex128(0).itemsize * M:
      dtype = numpy.complex128
    else:
      raise PyqmcDataError, "Invalid length of array data field: %s" % (data_rec_len,)
    return dtype

