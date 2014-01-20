#!/usr/bin/python
# $Id: pwqmc_meas.py,v 1.10 2010-09-30 01:45:42 wirawan Exp $
#
# pwqmc_meas.py
# Tools related to PWQMC-77 measurement dump (pwaf-*.ene) files
#
# Wirawan Purwanto
# Created: 20081022
#

# Python standard modules
import math
import os
import sys
import time
import weakref

#import tables # HDF5-based tables--nice, but not quite...
import h5py # HDF5 Python API
import numpy

from wpylib.timer import timer

class raw_meas_rec(object):
  """Raw measurement record from one of the "load" functions.
  This class holds the basic record, which contains the array of El and
  wtwlkr data.
  Generally this holds the measurement data from a population snapshot
  (or, "row").

  Note: here we only define El and wtwlkr as the member of the fields;
  but raw_meas_hdf5.load_meas_data also store `meta' field."""
  def __init__(self, El, wtwlkr):
    self.El = El
    self.wtwlkr = wtwlkr
    self.ndata = len(El)
  def __len__(self):
    return self.ndata
  # FIXME: add deepcopy() to return a deep copy [including the
  # metadata]


class raw_meas_view_blk(dict):
  """View of raw measurement records based on iblk/imeas hierarchy."""
  #def __init__(self):
  def add(self, iblk, imeas, rec):
    '''Adds a measurement record to the view.'''
    ## BUG: if existing iblk,imeas combination exists, then the block
    ## view is no longer valid!
    self.setdefault(iblk, {})
    self[iblk][imeas] = rec
    return rec
  # FIXME: add deepcopy() to return a deep copy


class raw_meas_view_flat(list):
  """View of raw measurement records in a flat 1-D array."""
  def add(self, rec):
    self.append(rec)
    return rec
  # FIXME: add deepcopy() to return a deep copy


# Predefined batch operations:
class batch_op(object):
  nop_global = lambda x: None  # a no-op for global level
  nop_global = staticmethod(nop_global)
  gather_global = lambda x: x  # a pass-through filter to get the intermediate list
  gather_global = staticmethod(gather_global)

  @staticmethod
  def max_field(fld):
    '''Returns Highest element in the fld field. fld is either El or wtwlkr.'''
    return (
      ( lambda rec: numpy.max(getattr(rec, fld)) ),
      ( lambda rec: numpy.max(rec) )
    )
  @staticmethod
  def min_field(fld):
    '''Returns Highest element in the fld field. fld is either El or wtwlkr.'''
    return (
      ( lambda rec: numpy.min(getattr(rec, fld)) ),
      ( lambda rec: numpy.min(rec) )
    )

#batch_op.nop_global = lambda x: None  # a no-op for global level


class raw_meas_data(object):
  """Object for holding raw measurement records returned by one of
  the "load" functions of the raw_meas_hdf5 object.
  This object holds an in-memory array of the data; updating the data will
  not update the HDF5 file from which the data has been taken.

  There are two ways of accessing the data: one is using iblk/imeas view,
  and another using flat view (no blocking)."""
  def __init__(self, blkview=True):
    """
    If blkview is False, then no hierarchical iblk/imeas view is possible.
    """
    if blkview:
      self.blk = raw_meas_view_blk()
    self.flat = raw_meas_view_flat()

  def add(self, iblk, imeas, El, wtwlkr):
    rec = raw_meas_rec(El, wtwlkr)
    if "blk" in dir(self):
      self.blk.add(iblk, imeas, rec)
    self.flat.add(rec)
    return rec

  def make_flat_section(self, Range):
    '''Returns a section of the raw_meas_data. Range is a slice object for the
    flat array.'''
    # FIXME: block view will be nice to add this later (not always possible)
    rslt = raw_meas_data(blkview=False)
    rslt.flat.extend(self.flat[Range]) # DIRTY HACK
    return rslt

  def __getitem__(self, idx):
    '''Accesses a given measurement record. For now, it supports the
    following ways of indexing:

    object[recno]  -- accesses the "recno"-th record in the flat view.
    object[iblk,:] -- accesses the "iblk"-th block in the blk view.
    object[iblk,imeas] -- accesses the [iblk][imeas] record in the blk view.

    Primitive slicing is supported, but the resulting record will lose its
    block hierarchy (FIXME):

    object[start:stop:step]  -- returns the corresponding slice of flat view,
          returned as another raw_meas_data object.
    '''
    if isinstance(idx, int):       # flat view
      return self.flat[idx]
    elif isinstance(idx, slice):   # slice of flat view
      return self.make_flat_section(idx)
    elif isinstance(idx, tuple):   # iblk,imeas view
      if "blk" not in dir(self):
        raise IndexError, \
          "Hierarchical iblk/imeas view is not possible with this viewer instance."
      (iblk,imeas) = idx
      if isinstance(imeas, slice): # disregarded
        # Instead of returning the lower-level object, we will
        # try to return a view object of same kind (raw_meas_data) here.
        #return self.blk[iblk]
        rslt = raw_meas_data(blkview=True)
        imeases = self.blk[iblk].keys()
        imeases.sort()
        for im in imeases:
          rec = self.blk[iblk][im]
          rslt.add(iblk, im, rec.El, rec.wtwlkr)
        return rslt
      elif isinstance(imeas, int):
        return self.blk[iblk][imeas]
    raise IndexError, "Invalid way of referring to a measurement record."

  def ndata(self):
    """Returns the number of records in the flat view."""
    return len(self.flat)

  # Two aliases:
  __len__ = ndata
  nrows = ndata

  def nblk(self):
    """Returns the number of blocks in the blocked view."""
    if "blk" not in dir(self):
      raise RuntimeError, \
        "Hierarchical iblk/imeas view is not possible with this viewer instance."
    return len(self.blk)

  def __iter__(self):
    """Returns iterator over the flat records, for simple iteration over
    the measurement records."""
    return self.flat.__iter__()

  def shape_stats(self):
    """Returns statistical information about the data shape.
    Returned as 5-tuple (or more, later):
    . number of records (rows)
    . average, stddev, min, and max of the number of measurement data in a row
    """
    popsizes = numpy.array([ len(row) for row in self ])
    return (len(self), popsizes.mean(),
            popsizes.std(),
            popsizes.min(), popsizes.max())

  def batch_operate(self, ranges, op, global_op = None):
    """Performs batch/collective operation on a range of records as
    given in ranges.
    Then we can do a final operation on the intermediate results.

    The per-record operator is given by op, and the global operator is
    given by global_op. Alternatively, op can be a 2-tuple containing
    the record- and global-level operators. """

    tmp_rslt = []
    if type(op) == tuple:
      global_op = op[1]
      op = op[0]
    if global_op == None: global_op = batch_op.nop_global
    for idx in ranges:
      tmp_rslt.append( op(self[idx]) )
    return global_op(tmp_rslt)

  def all_recs(self):
    """A shortcut to refer to all records in the flat view. For use with
    batch_operate."""
    return xrange(len(self.flat))

  # Static method for measurement data
  #@staticmethod
  rec_weighted_avg = \
      lambda rec: numpy.average(rec.El, weights=rec.wtwlkr, returned=True)

  @staticmethod
  def weighted_avg_op(fld, result = "array"):
    """Batch operator for calculating the weighted averages of the measurement
    records.
    The specific field to be weight-averaged must be specified in the "fld" argument
    (use "El" for local energy).
    Depending on the "result" setting, it can yield one of the following:
    - "array": a structured numpy array containing the weight-averages from each record
    - "wavg":  a tuple containing the global weighted average and the sum of all the
      weights (wtwlkr)
    - "wavg+err": like "wavg", but there is an additional field in the tuple, which is
      the error estimate of the mean.

    Example usage for typical weighted average over all measurement records
    (rows), supposing rec is a `raw_meas_data' object instance:

    rec.batch_operate(rec.all_recs(), rec.weighted_avg_op('El'))

    """
    if result == "array":
      global_op = \
          lambda arr: numpy.array(arr, dtype=[(fld, "=f8"), ("wtwlkr", "=f8")])
    elif result == "wavg":
      global_op = \
          lambda arr: numpy.average([ a[0] for a in arr ],
                                    weights=[ a[1] for a in arr ],
                                    returned=True)
    elif result == "wavg+err":
      global_op = \
          lambda arr: numpy.average([ a[0] for a in arr ],
                                    weights=[ a[1] for a in arr ],
                                    returned=True) + \
                      (numpy.std([ a[0] for a in arr ]) / math.sqrt(len(arr) - 1),)
      # Using sqrt(len(arr)-1) to compensate the biased estimate in numpy.std ^^^
    else:
      raise ValueError, "Invalid result type for weighted_avg_op: " + result

    return (
      lambda rec: numpy.average(getattr(rec,fld), weights=rec.wtwlkr, returned=True),
      global_op
    )

  @staticmethod
  def cap_Elocal_op(El_caps, verbose = False):
    """Batch operator to cap the local energies in the measurement records.
    Yields None."""
    meta = {"El_caps" : El_caps, "verbose" : verbose}
    return (
      (lambda rec: cap_Elocal(rec.El, **meta)),
      batch_op.nop_global
    )



class raw_meas_hdf5(object):
  """Creates a new accessor for a `raw measurement' data structure.
  This object is not to be created directly by users; instead, use
  meas_hdf5.raw_open() object method to create the accessor.

  NOTE: the accessor object's validity depends on the life of the parent
  meas_hdf5 object."""
  def __init__(self, raw_group, parent):
    """Creates a new raw measurement accessor.
    This routine is to be called by the meas_hdf5.raw_open() object method.
    """
    # The reference to parent *MUST* be a weak reference or else it would
    # cause circular reference that prevents the parent object from being
    # deleted properly.
    self.parent = weakref.ref(parent)
    self.raw_group = raw_group
    self.E_l = raw_group['E_l']
    self.wtwlkr = raw_group['wtwlkr']
    self.proc = raw_group['proc']
    # Measurement metadata:
    self.meta = raw_group['meta']

  def check_alive(self):
    '''Checks if the parent object is still alive and valid.'''
    if self.parent == None or self.parent() is None:
      raise ReferenceError, \
          "The parent meas_hdf5 object has died or been closed"

  def job(self):
    '''Returns the HDF group of the associated QMC job (which is the parent
    directory of this raw group).'''
    job_group = self.raw_group.name
    job_group = job_group[0:job_group.rfind("/")]
    return self.parent().dbh5[job_group]

  def append_row(self, iblk, imeas, beta, \
                 stage = 2, select_row = False, \
                 reserve_size = None):
    '''Appends a new row on the measurement record.
    The select_row option makes the newly created row as the current row.
    This method primarily records only the new row's metadata into the
    metadata table.
    The `reserve_size' option can be used to make room for at least
    `reserve_size' number of elements in the measurement data storage.
    Returns the integer index of the new measurement row.'''
    self.check_alive()
    g = self.raw_group

    do_extend = True
    shape = self.E_l.shape
    new_row_size = shape[0]+1
    new_col_size = shape[1]
    if shape[0] == 1:
      # Special condition: if there is only 1 row and its count field is
      # negative, than it is a dummy empty element.
      # In this case we should not extend the shape unless the reserve
      # size requires that.
      meta_old = self.meta[0]
      if meta_old[3] < 0:
        do_extend = False
        new_row_size = shape[0]

    # If we reserve a larger column size, we always need to extend:
    if reserve_size > new_col_size:
      new_col_size = int(reserve_size)
      do_extend = True

    if do_extend:
      shape = self.extend_shape((new_row_size, new_col_size))
    new_row = shape[0]-1

    if reserve_size != None and reserve_size > 0:
      # This is a trick to fill the chunked sectors to at least
      # reserve_size elements so that it lies contiguous (hopefully)
      # in the HDF5 file:
      # FIXME: This can get VERY LARGE if reserve_size is huge.
      # We should use some kind of iterator in that case.
      #print "reserving ", reserve_size
      Zeros = numpy.zeros((1, reserve_size), dtype="=i4")
      self.E_l[new_row, 0:reserve_size] = Zeros.astype("=c16")
      self.wtwlkr[new_row, 0:reserve_size] = Zeros.astype("=f8")
      self.proc[new_row, 0:reserve_size] = Zeros

    # Initialize the new metadata:
    meta_new = numpy.array((iblk, imeas, beta, 0, stage), \
                           dtype=self.parent().raw_meta_dtype)
    self.meta[new_row] = meta_new

    if select_row: self.row = new_row
    return new_row

  def extend_shape(self, new_shape):
    '''Extends the shape of the data structure.'''
    self.check_alive()
    self.E_l.id.extend(new_shape)
    self.wtwlkr.id.extend(new_shape)
    self.proc.id.extend(new_shape)
    self.meta.id.extend(new_shape[0:1])
    return new_shape

  def shape(self):
    '''Returns the current shape of the data structure.'''
    self.check_alive()
    return self.E_l.shape

  def row_count(self):
    self.check_alive()
    return self.E_l.shape[0]

  def select_row(self, row):
    '''Selects a particular measurement row to work on (to fetch data, etc.)'''
    self.check_alive()
    row_count = self.row_count()
    if row < 0 or row >= row_count:
      raise IndexError, 'Invalid row index for raw measurement data'
    else:
      self.row = row

  def build_row_cache(self, *criteria):
    '''Builds a row cache for fast seek later. It is using simple dict to
    map criteria to the flat indices. If a key is not unique, then only
    the first occurence of each key will be used.'''
    self.check_alive()
    fld_index = {}
    i = 0
    for key in self.parent().raw_meta_dtype.names:
      fld_index[key] = i
      i += 1
    flds = [ fld_index[key] for key in criteria ]
    #print flds

    row = 0
    cache = {}
    for m in self.meta:
      #print "key", ",".join([ str(m[key]) for key in flds ]), " = ", row
      cache[ ",".join([ str(m[key]) for key in flds ]) ] = row
      row += 1
    return cache

  def find_row(self, _create_new = None, _select_row = False, **criteria):
    '''Finds a row matching certain criteria. The criteria is given as
    "kw1=value1, kw2=value2, ..." variable arguments.

    The optional _create_new argument is there for specifying additional
    parameters for creating a new row, in case no row matches the given
    search criteria.
    '''
    self.check_alive()
    meta = self.parent().raw_meta_blank()
    #print "find_row: looking for", criteria
    # First, make a one-level copy to save the original hash:
    criteria = criteria.copy()
    # Extract the first criterion for quick comparison; do the rest in
    # iterative way:
    for (crit1, val1) in criteria.iteritems():
      del criteria[crit1]
      break
    row = 0
    for m in self.meta:
      #print " ", m
      meta[0] = m
      if meta[crit1] == val1:
        Matched = True
        for (c, v) in criteria.iteritems():
          if meta[c] != v:
            Matched = False
            break
        if Matched:
          #print "find_row: matched @ row", row
          if _select_row: self.select_row(row)
          return row
      row += 1

    # If we arrive here, then we haven't found the record yet. See if _create_new
    # is set:
    if _create_new != None:   # then it *must* be a hash
      criteria[crit1] = val1  # re-add the first criterion
      criteria.update(_create_new) # update with more values for creating new stuff
      row = self.append_row(select_row = _select_row, **criteria)
      #print "find_row: creating a new row ", row, " with ", criteria
      return row
    else:
      return None # in case there is no match

  def row_append_meas(self, E_l, wtwlkr, proc):
    '''Appends measurement data into the current row. All input arrays must have
    the same dimensionality.'''
    # TODO FIXME: Remove the staging arrays if E_l and all the other arrays are
    # already the correct numpy array with same dtype.
    self.check_alive()
    E_l_ = numpy.array(E_l, dtype=self.E_l.dtype)

    if getattr(wtwlkr, "__iter__", False):
      wtwlkr_ = numpy.array(wtwlkr, dtype=self.wtwlkr.dtype)
      if E_l_.shape != wtwlkr_.shape:
        raise TypeError, "wtwlkr array has different dimensionality from E_l"
    else:
      # Assume scalar; use broadcasting
      wtwlkr_ = numpy.empty(E_l_.shape, dtype=self.wtwlkr.dtype)
      wtwlkr_[:] = wtwlkr

    if getattr(proc, "__iter__", False):
      proc_ = numpy.array(proc, dtype=self.proc.dtype)
      if E_l_.shape != proc_.shape:
        raise TypeError, "proc array has different dimensionality from E_l"
    else:
      # Assume scalar
      proc_ = numpy.empty(E_l_.shape, dtype=self.proc.dtype)
      proc_[:] = proc

    # Expand the row if necessary:
    col_size = self.row_ndata()
    col_free = self.E_l.shape[1] - col_size
    col_needed = E_l_.shape[0] - col_free
    if col_needed > 0:
      new_size = (self.E_l.shape[0], self.E_l.shape[1] + col_needed)
      self.extend_shape(new_size)

    row = self.row
    new_begin = col_size
    new_end = col_size + E_l_.shape[0]
    self.E_l[row, new_begin:new_end] = E_l_
    self.wtwlkr[row, new_begin:new_end] = wtwlkr_
    self.proc[row, new_begin:new_end] = proc_
    # Update the 'count' metadata: cannot use simple assignment for this
    new_meta = numpy.array(self.meta[row], dtype=self.parent().raw_meta_dtype)
    new_meta['count'] += E_l_.shape[0]
    self.meta[row] = new_meta

  def row_ndata(self, row=None):
    '''Queries the number of data points existing in a row.'''
    self.check_alive()
    if row == None: row = self.row
    return self.meta[row][3]

  def load_meas(self, blocks, indices, stage = 2):
    '''A master routine to load measurement across many blocks/indices
    into memory (compatible with early "slices" format).
    Returns the snapshots of (1) wtwlkr and (2) local energy as a tuple.

    Be aware that if the measurement files are large, then the data
    structure in memory is also very large!

    DEPRECATED ROUTINE. Please use load_meas_data instead!'''

    self.check_alive()
    wtwlkr_slices = {}
    El_slices = {}
    cache = self.build_row_cache("stage", "iblk", "imeas")

    for blk in blocks:
      wtwlkr_slices.setdefault(blk, {})
      El_slices.setdefault(blk, {})
      for idx in indices:
        key = str(stage) + "," + str(blk) + "," + str(idx)
        self.select_row(cache[key])
        ndata = self.row_ndata()
        El_slices[blk][idx] = numpy.real(self.E_l[self.row, 0:ndata])
        wtwlkr_slices[blk][idx] = self.wtwlkr[self.row, 0:ndata]
      #print rslt
    return (wtwlkr_slices, El_slices)

  def load_meas2(self, stage = None):
    '''Alternate routine to load ALL measurements for a particular stage
    without regard to the block/index selection.

    DEPRECATED ROUTINE. Please use load_meas_data instead!'''
    self.check_alive()
    wtwlkr_slices = {}
    El_slices = {}
    meas_recs = []
    for row in xrange(0, self.row_count()):
      m = self.meta[row]
      #print m
      iblk = m[0]
      imeas = m[1]
      ndata = m[3]
      if ndata == 0: continue # don't add it in
      if stage != None and m[4] != stage: continue
      wtwlkr_slices.setdefault(iblk, {})
      El_slices.setdefault(iblk, {})
      El_slices[iblk][imeas] = numpy.real(self.E_l[row, 0:ndata])
      wtwlkr_slices[iblk][imeas] = self.wtwlkr[row, 0:ndata]
      meas_recs += [
        (
         El_slices[iblk][imeas],
         wtwlkr_slices[iblk][imeas],
         m,
        )
      ]
    return (wtwlkr_slices, El_slices, meas_recs)

  def load_meas_data(self, stage=None, match=None, append_to_raw=None,
                     value_processor=None,
                     keep_El_imag=True, keep_phasefac=True):
    """Loads measurement records for a given stage (if `stage' is specified) or
    those matching given record (if `match' is specified).

    The `match' function must accept a dict argument, which contains metadata
    of the row under consideration.
    It must return True if the row is to be accepted; False otherwise.
    Note that `stage' and `match' should not be specified simultaneously.
    Returns a raw_meas_data object for convenience.

    The `append_to_raw' argument, if given, must be an existing raw_meas_data
    object to which we want to append the dataset.
    This option is useful, for example, for chaining the raw data from various
    sequential parts of a QMC run.

    The `value_processor', if given, must be a python function with the
    following prototype:
        value_processor(<value_dict>, <metadata_dict>)
    The values are given as <value_dict> with three keys:
        'E_l', 'wtwlkr', 'proc'.

    The arguments `keep_El_imag' and `keep_phasefac' have similar meanings
    as in the conversion routines.
    By default these are True, which means that the imaginary part of E_l and
    the complex phase of wtwlkr are included, unless if they have already been
    discarded when converting from the *.ene files (see notes in
    convert_meas_to_hdf5_v2 routine).
    """
    if append_to_raw != None:
      rslt = append_to_raw
    else:
      rslt = raw_meas_data()
    # This way does not give memory advantage:
    #El = numpy.empty((self.E_l.shape[1],), dtype="=f8")
    rows_selected = []
    metas = []
    for row in xrange(0, self.row_count()):
      m = dict(zip(self.parent().raw_meta_dtype.names, self.meta[row]))
      #print m
      if m['count'] == 0: continue # don't add record w/ zero length
      m['row'] = row
      if stage != None and m['stage'] != stage: continue
      if match != None and not match(m): continue
      rows_selected.append(row)
      metas.append(m)

    for (row,m) in zip(rows_selected, metas):
      iblk = m['iblk']
      imeas = m['imeas']
      ndata = m['count']
      #El[0:ndata] = self.E_l[row, 0:ndata].real

      E_l = self.E_l[row, 0:ndata]
      wtwlkr = self.wtwlkr[row, 0:ndata]

      if value_processor != None:
        proc = self.proc[row, 0:ndata]
        meta = m.copy()
        meta['row'] = row
        data = { 'E_l': E_l, 'wtwlkr': wtwlkr, 'proc': proc }
        value_processor(data, meta)

      if keep_El_imag:
        val_El = E_l
      else:
        El=numpy.copy(E_l.real)
      if keep_phasefac:
        val_wtwlkr = wtwlkr
      else:
        val_wtwlkr = numpy.abs(wtwlkr)
      rec = rslt.add(iblk, imeas, \
                     El=val_El, \
                     wtwlkr=val_wtwlkr)
      # Also store the metadata here:
      rec.meta = m

    return rslt

  def get_raw_object(self, src_raw):
    # No support for derived class yet :-(
    if isinstance(src_raw, meas_hdf5):
      src_raw = src_raw.raw
    elif isinstance(src_raw, h5py.highlevel.Group):
      src_raw = src_raw['raw']
    elif isinstance(src_raw, raw_meas_hdf5):
      pass
      #src_raw = src_raw.raw_group --nope
    else:
    #elif type(src_raw) != raw_meas_hdf5:
      raise ValueError, \
          "Don't know how to handle object of type " + str(type(src_raw))
    return src_raw

  def compare_metadata(self, src_raw, verbose=0):
    '''Compares metadata from another raw_meas_hdf5 object to see if it
    matches the signature of this object.'''

    mismatch = 0
    src_raw = self.get_raw_object(src_raw)
    thisjob = self.job()
    thatjob = src_raw.job()
    thismeta = thisjob.attrs
    thatmeta = thatjob.attrs
    for key in ['System', 'H0', 'deltau', 'Evar', 'units', 'kptstr',
                'GeomCode', 'wsvol']:
      if thismeta[key] != thatmeta[key]:
        mismatch += 1
        if verbose:
          if mismatch == 1:
            print "Mismatch(es) found in raw_meas_hdf5 metadata (this <-> other):"
          print "*", key, ":", thismeta[key], "<=>", thatmeta[key]
    return mismatch

  def append_raw_records(self, src, warn_meta_mismatch = None, debug=0):
    '''Appends raw measurement records from another database's raw records.
    The src must be a raw_meas_hdf5 object.'''

    src = self.get_raw_object(src)

    if warn_meta_mismatch == None:
      warn_meta_mismatch = sys.stdin.isatty()

    if debug >= 1:
      print "Appending raw records from file ", src.parent().filename(), \
            "to", self.parent().filename()

    Cmp = self.compare_metadata(src, verbose=1)
    if not Cmp:
      pass # OK, can append
    elif warn_meta_mismatch:
      print "This database = ", self.parent().filename()
      print "That database = ", src.parent().filename()
      print "append_raw_meas_data: Do you still want to continue appending records? [y/n]",
      ans = sys.stdin.readline().strip().lower()
      if not ans.startswith("y"):
        return False
    else:
      raise ValueError, "Metadata mismatch"

    #beta0 = self.meta['beta'][self.row_count()-1]
    meta_names = self.parent().raw_meta_dtype.names

    jobmeta = self.job().attrs
    itv_Em = jobmeta['itv_Em']
    nblkstep = jobmeta['nblkstep']
    last_dest_meta = self.meta[self.row_count()-1]
    dest_m = dict(zip(meta_names, last_dest_meta))
    beta0 = dest_m['beta']
    iblk = dest_m['iblk']
    imeas = dest_m['imeas']
    for src_row in xrange(src.row_count()):
      src.select_row(src_row)
      src_meta = src.meta[src_row]
      m = dict(zip(meta_names, src_meta))  # -- can't update new meta
      # FIXME: must be able to start at an existing dest row (i.e. "overwrite
      # mode"), rather than always appending at the end.
      # FIXME: must be able to delete rows if final(ndata) < nrows (dest original value)
      # FIXME: for now we assume the appended recs are from measurement phase
      # thus we continue the blocking param rather than using the existing ones
      # FIXME: still can't do incommensurate measurement interval
      src_ndata = m['count']

      imeas += itv_Em
      if imeas > nblkstep:
        iblk += 1
        imeas -= nblkstep

      if debug >= 100:
        print "append_row #", src_row, ":", \
            (iblk, imeas, m['beta'] + beta0, m['stage'], src_ndata)
      self.append_row(iblk, imeas, \
                      beta=m['beta'] + beta0, \
                      stage=m['stage'], \
                      select_row=True, reserve_size=src_ndata)

      self.row_append_meas(src.E_l[src_row, 0:src_ndata],
                           src.wtwlkr[src_row, 0:src_ndata],
                           src.proc[src_row, 0:src_ndata])

    self.parent().flush()



class meas_hdf5(object):
  """The main object for storing PWQMC-77 measurement results in HDF5 format.
  """
  # These static fields are put here instead of in the raw_meas_hdf5
  # class, thus permitting customization of the parameters as necessary:
  raw_meta_dtype = numpy.dtype([
    ('iblk', '=i4'),
    ('imeas', '=i4'),
    ('beta', '=f8'),
    ('count', '=i4'),
    ('stage', '=i2'),
  ])
  # Default datatypes for raw structure:
  raw_dtypes = numpy.dtype([
    ("E_l", "=c16"),
    ("wtwlkr", "=f8"),
    ("proc", "=i4"),
    ("meta", raw_meta_dtype),
  ])
  # The default number of data is -1 (<0) as a special marker for the first
  # blank row.
  # stages: 0 = eqlb, 1 = growth, 2 = measurement
  meta_blank0 = (0, 0, 0.0, -1, 2)
  default_raw_group = '/job0/raw'
  # default chunking: 4 in time direction, 32 in wlkr direction
  default_raw_chunks = [4, 32]
  # default header metadata:
  header_meta = {
    "Creator": "pyqmc.results.pwqmc_meas",
    "FormatType": "pwqmc_meas",
    "FormatVersion": 0,
    "Devel_CVS_ID": "$Id: pwqmc_meas.py,v 1.10 2010-09-30 01:45:42 wirawan Exp $",
  }
  # FIXME:
  # - delete Devel_CVS_ID
  # - add date/time of creation
  # - ...what else?

  def __init__(self, fname=None, mode="a", create_raw=True):
    if fname: self.open(fname, mode=mode, create_raw=create_raw)

  def __del__(self):
    #if hasattr(self, "dbh5")
    # We must delete circular reference to itself via raw ptr:
    #if hasattr(self, "raw"):
    #  self.raw.parent = None
    #  self.raw = None
    self.close()

  def open(self, fname, mode="a", create_raw=True):
    '''Opens or creates a HDF5 database containing the PWQMC measurement dump.
    This is a more efficient way of storing the pwaf-*. raw data.'''
    self.close()
    self.dbh5 = h5py.File(fname, mode)
    if "FormatType" not in self.dbh5.attrs:
      self.create_new(create_raw)
    try:
      self.raw_open(self.default_raw_group)
    except:
      pass

  def close(self, debug=1):
    '''Closes a previously opened HDF5 measurement database.'''
    if hasattr(self, "dbh5"):
      if debug: print "closing hdf5 database " + self.filename()
      self.dbh5.close()
      delattr(self, "dbh5")
    # We must also close circular reference to itself via "raw" reference:
    if hasattr(self, "raw"):
      self.raw.parent = None
      delattr(self, "raw")

  def flush(self):
    self.dbh5.flush()

  def filename(self):
    """Returns the filename of the underlying HDF5 database.
    This is given to prevent confusion in the API to get the filename:
    version 1.0 and 1.1 use 'name' attribute but newer versions use
    'filename' attribute."""
    if h5py.version.version_tuple >= (1,2,0):
      return self.dbh5.filename
    else:
      return self.dbh5.name

  def create_new(self, create_raw=True):
    '''Creates a new HDF5 structure for a measurement result.'''
    db = self.dbh5
    for (k,v) in self.header_meta.iteritems():
      db.attrs[k] = v
    #db.attrs["Creator"] = "pyqmc.results.pwqmc_meas"
    db.attrs["CreateTime"] = \
        time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
    #db.attrs["FormatVersion"] = 0
    #db.attrs["FormatType"] = "pwqmc_meas"
    if create_raw:
      self.create_group_path(self.default_raw_group)
      self.raw_create(self.default_raw_group)
    db.flush()

  def create_group_path(self, path):
    '''Creates HDF5 group and subgroups that form the path.
    The path must be an absolute path.'''
    if not path.startswith('/'):
      raise ValueError, 'The "path" argument must be an absolute path!'
    paths = path.split('/')
    # ^ don't include the first part, which must be a slash.
    # ^ don't include the last (group) part
    for lpath in xrange(2, len(paths)+1):
      g = self.dbh5.require_group("/".join(paths[0:lpath]))
    return g

  def raw_create(self, path, dtypes={}):
    '''Creates a new raw measurement section at absolute group location `path'.'''
    g = self.create_group_path(path)
    '''
    Notes for arrays E_l, wtwlkr, and proc:
    * these all must have identical dimensions
    * index #1 `meas_idx': global measurement index (0, 1, 2, 3, ...)
    * index #2 `wlkr_idx': walker index (0, 1, 2, 3, ..., meta[meas_idx].count)

    The meta array runs with meas_idx only.
    '''

    def get_dtype(fld):
      try:
        return dtypes[fld]
      except:
        return self.raw_dtypes[fld]

    Chunk = tuple(self.default_raw_chunks)
    g.create_dataset("E_l", shape=(1, 1), dtype=get_dtype("E_l"), \
        maxshape=(None,None), chunks=Chunk)
    g.create_dataset("wtwlkr", shape=(1, 1), dtype=get_dtype("wtwlkr"), \
        maxshape=(None,None), chunks=Chunk)
    g.create_dataset("proc", shape=(1, 1), dtype=get_dtype("proc"), \
        maxshape=(None,None), chunks=Chunk)
    # Measurement markers (blk/imeas/ compound):
    meta_blank = self.raw_meta_blank()[0]
    g.create_dataset("meta", shape=(1,), dtype=get_dtype("meta"), \
        maxshape=(None,), chunks=Chunk[0:1])
    g['meta'][0] = meta_blank
    return g

  def raw_meta_blank(self):
    return numpy.array([self.meta_blank0], dtype=self.raw_meta_dtype)

  def raw_open(self, path = None):
    '''Opens an existing raw measurement section for work.
    This provides convenient shortcuts (E_l, wtwlkr, etc.) at the object
    instance level.'''
    if path == None: path = self.default_raw_group
    g = self.dbh5[path]
    raw = raw_meas_hdf5(g, self)
    self.raw = raw  # last opened raw group
    # (WARNING: circular reference is in raw_meas_hdf5 object; to avoid
    # circular reference woes, we make the meas_hdf5 reference in the
    #  raw object "weak".)
    return raw


# end class meas_hdf5



class meas_text(object):
  '''meas_text is a class for reading the measurement dump produced by
  PWQMC-77 code (i.e. a single pwaf-*.meas or pwaf-*.ene file per MPI
  process).
  The raw dump to be read is in plain text format.

  Additional fields possessed by this object:
  * proc_rank = MPI rank (read from the dump file)
  * start_date, start_time = start date and time of the run (read from the
    dump file)
  * lastpos = last position of the file pointer. Will be changed whenever
    (1) a file is just opened, (2) a seek to a section, or (3) a read to
    a measurement section.
  * ndata = number of measurement points read by the last read_section call.
  '''

  def __init__(self, fname = None):
    if fname: self.open(fname)
  def __del__(self):
    self.close()

  # Opening and closing input files:
  def open(self, fname):
    '''Opens a given measurement file for reading and analysis.
    '''
    self.close()
    self.file = open(fname, "r")
    #self.fname = fname  -- use self.file.name instead
    # Requires a proper header:
    F = self.file
    self.check_file_header(F)

    while True:
      lastpos = F.tell()
      L = F.readline().strip()
      if L.startswith("#"): # skip comments
        pass
      elif F == "":
        raise ValueError, "File `" + fname + "': no measurement found"
      elif L.startswith("Process rank"):
        self.proc_rank = int(L.split()[2])
      elif L.startswith("Start time, date:"):
        (self.start_date, self.start_time) = L.split()[3:5]
      elif L.startswith("MEA ") or L.startswith("MEAS ") or \
           L.startswith("GRTH ") or L.startswith("EQLB "):
        F.seek(lastpos)
        break
    self.lastpos = lastpos

  def close(self):
    '''Closes a previously opened measurement file.'''
    if hasattr(self, "file"):
      self.file.close()
      delattr(self, "file")
      #delattr(self, "lastpos") #!!!

  def check_file_header(self, F):
    L = F.readline().strip()
    # Support old dataformat for now (generated by hacked v17)
    if not L.startswith("PWQMC measurement data: total energies") and \
       not (L == "PWQMC measurement data"):
      raise ValueError, "File `" + F.name + "': Invalid file type"

  def rewind(self):
    self.file.seek(0)
    self.lastpos = self.file.tell()

  def seek_section(self, section, block, index):
    '''Seeks forward a named section with a particular block and index.'''
    sect = section + " "
    F = self.file
    line = F.readline()
    #for line in self.file:
    while line != "":
      # This kind of seek should be very fast:
      if line.lstrip().startswith(sect):
        flds = line.split()
        if int(flds[1]) == block and int(flds[2]) == index:
          self.lastpos = F.tell()
          return True
      line = F.readline()
    return False

  def seek_any_section(self):
    '''Seeks forward to the next named section (any kind).
    Returns a tuple of 4 elements: (stage, block, index, ndata).
    Stage is an integer of value 0 (equilibration), 1 (growth), or 2
    (measurement). Block is the block number (0, 1, 2, ...).
    Index indicates the step counter (cf. variable 'istp') in the main
    random-walk loop at this particular block.
    Ndata is the number of data points in this measurement block (per
    MPI process).
    On the earliest (obsolete) dataset format, Ndata field is nonexistent,
    and it is set to None.
    This function returns None if there is no more measurement section
    detected.
    '''
    F = self.file
    line = F.readline()
    #for line in self.file:
    while line != "":
      # This kind of seek should be very fast:
      L = line.lstrip()
      if L.startswith('MEAS ') or L.startswith('MEA ') or \
         L.startswith('EQLB ') or L.startswith('GRTH '):
        flds = line.split()
        if L.startswith('MEAS ') or L.startswith('MEA '):
          # MEAS is the recommended keyword; MEA is old keyword
          stage = 2
        elif L.startswith('GRTH '):
          stage = 1
        else:
          stage = 0
        block = int(flds[1])
        index = int(flds[2])
        if len(flds) >= 4:
          ndata = int(flds[3])
        else:
          ndata = None
        self.lastpos = F.tell()
        return (stage, block, index, ndata)
      line = F.readline()
    return None

  def read_section(self, section=None, block=1, index=1, \
    fname=None, keep_El_imag=False, keep_phasefac=False, H0=0, out=None):
    '''Reads energy measurement, one section at a time from a pwaf-*.ene file.

    On entry:
    * If section is set to None, then it begins reading the measurement section
      immediately.
      Otherwise, section, block, and index arguments must be specified, and
      the a forward seek to the specified section will be performed beforehand.
    * if keep_El_imag is True, then the imaginary part of E_l is preserved
      (thus E_l must be a complex array).
    * if keep_phasefac is True, then the phase factor is also read (this is only
      for free projection runs).
    * H0 is the constant part of the energy.
    * out, is given, must be a structured numpy array that will contain the
      result. The dtype must contain:
      - ('wtwlkr',float)   (complex if keep_phasefac==True)
      - ('E_l',float)      (complex if keep_El_imag==True)

    On exit:
    * The file pointer is positioned to read the immediate next section.
    * Returns a 2-component tuple containing (1) wtwlkr and (2) local energy
      data. If out is given, then tuple components point to the corresponding
      fields in out.

    Valid section names are "EQLB", "GRTH", or "MEAS" (according to the section
    strings dumped by PWQMC code).
    '''
    from wpylib.math import complex_polar
    if fname != None: self.open(fname)
    if section != None and not self.seek_section(section, block, index):
      raise ValueError, \
        "Cannot seek section " + section + " block " \
        + str(block) + " index " + str(index)

    if out == None:
      # Legacy data structure: using lists
      # This is deprecated; please use the numpy array buffer, which
      # is better!
      wtwlkr = []
      E_l = []
    else:
      wtwlkr = out['wtwlkr']
      E_l = out['E_l']

    count = 0
    F = self.file
    lastpos = F.tell()
    line = F.readline()
    self.next_section = None
    while (line != ""): # for line in self.file:
      #print "readline:%d: %s" % (lastpos, line)
      flds = line.split()
      if not flds[0].isdigit():
        # if next section is encountered, undo the last line read:
        self.file.seek(lastpos)
        self.next_section = flds
        break

      if keep_phasefac:
        w = complex_polar(float(flds[1]), float(flds[4]))
      else:
        w = float(flds[1])
      if keep_El_imag:
        E = complex(float(flds[2]), float(flds[3])) + H0
      else:
        E = float(flds[2]) + H0

      if out == None:
        wtwlkr.append(w)
        E_l.append(E)
      else:
        wtwlkr[count] = w
        E_l[count] = E

      count += 1
      lastpos = F.tell()
      line = F.readline()

    self.lastpos = lastpos
    self.ndata = count
    rslt = (wtwlkr, E_l)
    return rslt

  def make_numpy_buffer(self, size, keep_El_imag=False, keep_phasefac=False):
    '''Creates a numpy structured array for read_section routine's `out'
    argument.'''
    if keep_El_imag:
      El_type = complex
    else:
      El_type = float

    if keep_phasefac:
      wt_type = complex
    else:
      wt_type = float

    if isinstance(size, tuple):
      sz = size
    else:
      sz = (size,)

    return numpy.empty(sz, dtype=[('wtwlkr',wt_type), ('E_l',El_type), ('proc','=i4')])



def append_array(arr, dest, idx):
  '''Private routine to append a list at the end of a one-dimensional array.'''
  if idx not in dest:
    dest[idx] = numpy.array(arr)
    #print "created, idx = ", idx, "size = ", dest[idx].shape
  else:
    dest[idx] = numpy.concatenate((dest[idx], arr))
    #print "appended, idx = ", idx, "size = ", dest[idx].shape


def load_measurements(blocks, indices, H0 = 0, files = None):
  '''A master routine to load measurement across many pwaf-*.meas files into
  memory.
  Returns the snapshots of (1) wtwlkr and (2) local energy as a tuple.

  Be aware that if the measurement files are large, then the data structure in
  memory is also very large!

  DEPRECATED ROUTINE. It is preferrable to convert the *.ene files to HDF5 format
  then work with the HDF5 database. The analysis tools are more extensively
  developed for that format!'''

  files = meas_glob(files)
  wtwlkr_slices = {}
  El_slices = {}
  mea_f = meas_text()
  for fname in files:
    #if f > "pwaf-00001.meas": break
    print "Reading", fname
    mea_f.open(fname)
    for blk in blocks:
      #if not (blk in wtwlkr_slices): wtwlkr_slices[blk] = {}
      #if not (blk in El_slices): El_slices[blk] = {}
      wtwlkr_slices.setdefault(blk, {})
      El_slices.setdefault(blk, {})
      for idx in indices:
        (wtwlkr1, El1) = mea_f.read_section("MEA", blk, idx, H0=H0)
        append_array(wtwlkr1, wtwlkr_slices[blk], idx)
        append_array(El1, El_slices[blk], idx)
    mea_f.close()
    #print rslt
  #wtwlkr = numpy.array(w)
  #E_l = numpy.array(E)
  return (wtwlkr_slices, El_slices)


def convert_meas_to_hdf5(output, H0=0, files=None, **opts):
  '''A master routine to convert measurement across many pwaf-*.meas files to
  the new hdf5 format.

  Additional options:
  . betablk: imaginary time length of a block
  . deltau: imaginary time step
  . optimize: finds out the number of measurements in all files
    and use that information to optimize the storage use and layout.
    This can be a very long operation if the files are large or many.
  . nmeas_counts: provides the number of measurement points in each record.
    This is an alternate way to optimize the file layout and size.
    Usually this is not necessary as it can be provided automatically via
    the optimize option.
  . match: optional subroutine to filter in (or out) the measurement blocks
    to be included (or excluded) from the measurement phase.
    The `match' function must accept a hash argument, which contains metadata
    of the row under consideration.
    It must return True if the row is to be accepted; False otherwise.

  Advanced options:
  . meas_text_class
  . meas_hdf5_class
    Used to define custom classes for creating measurement reader
    (meas_text_class) and archiver (meas_hdf5_class).
    Know what you are doing if you specify this.
    They must be API-compatible to the meas_text and meas_hdf5 defined
    in this module.

  NOTE: This routine is deprecated.
  Please consider using the more efficient convert_meas_to_hdf5_v2 routine.
  '''
  global hdf5_conv_last
  from pyqmc.stats.avg import avg

  files = meas_glob(files)

  T = timer()
  tm_meas_count = avg()

  meas_hdf5_class = opts.get("meas_hdf5_class", meas_hdf5)
  meas_text_class = opts.get("meas_text_class", meas_text)
  optDebug = opts.get("debug", 0)
  optTime = opts.get("time", 0)
  if opts.get('optimize', False):
    T.start()
    nmeas_counts = get_meas_count(files)
    tm_meas_count += T.stop()
  else:
    nmeas_counts = opts.get("nmeas_counts", None)

  mea_f = meas_text_class()
  rslt = meas_hdf5_class()
  # Must set the chunking parameter here if required:
  # Set chunk size to ~2% of the avg measurement count.
  if nmeas_counts != None:
    chunk_size = max(rslt.default_raw_chunks[1], \
                     int(numpy.mean(nmeas_counts.values())/50))
    rslt.default_raw_chunks = [1, chunk_size]
    if optDebug: print "New chunk size = ", rslt.default_raw_chunks

  rslt.open(output)
  raw_hdf5 = rslt.raw
  cache = raw_hdf5.build_row_cache('stage', 'iblk', 'imeas')
  # provides a handle to the last created db just in case the
  # conversion fails in the middle:
  hdf5_conv_last = rslt
  # Reads off optional parameters
  betablk = opts.get('betablk', 0.5)
  deltau = opts.get('deltau', 0.01)
  match = opts.get('match', None)

  tm_readsect = avg()
  tm_find = avg()
  tm_append = avg()
  tm_seek = avg()
  for fname in files:
    #if f > "pwaf-00001.meas": break
    #print "Reading", fname
    mea_f.open(fname)
    if optDebug:
      print "File ", fname, "rank = ", mea_f.proc_rank, \
            "start @ ", mea_f.start_date, mea_f.start_time
    #print "<" + line + ">"

    section = mea_f.seek_any_section()
    while section:
      (stage, iblk, imeas, ndata) = section
      if match:
        if not match({'rank': mea_f.proc_rank,
                      'stage': stage, 'iblk': iblk, 'imeas': imeas,
                      'count': ndata }):
          T.start()
          section = mea_f.seek_any_section()
          tm_seek += T.stop()
          continue

      #print "read_sect"
      T.start()
      meas = mea_f.read_section(section=None, H0=H0, keep_El_imag=True)
      tm_readsect += T.stop()
      #print "stage ", next_section[0], " blk ", next_section[1], " idx ", next_section[2]
      #print rslt[0]
      #print rslt[1]
      # The row MUST be created regardless whether this particular measurement
      # record is empty or not.
      # FIXME for beta below: this is temporary, only for commensurate
      # measurement interval:
      beta = (iblk - 1) * betablk + imeas * deltau
      T.start()
      '''# THIS WAY IS INEFFICIENT:
      row = raw_hdf5.find_row(iblk=iblk, imeas=imeas, stage=stage, \
                              _create_new = {'beta': beta}, \
                              _select_row = True)
      '''
      # For now we manage the cache ourselves. Maybe later we can use a
      # different means.
      key = ",".join( (str(stage), str(iblk), str(imeas)) )
      if key not in cache:
        extra_opts = {}
        if nmeas_counts != None:
          if key in nmeas_counts:
            extra_opts["reserve_size"] = nmeas_counts[key]
        row = raw_hdf5.append_row(iblk=iblk, imeas=imeas, beta=beta, \
                                  stage=stage, select_row=True, \
                                  **extra_opts)
        cache[key] = row
      else:
        raw_hdf5.select_row(cache[key])

      tm_find += T.stop()
      if len(meas[0]) > 0:
        #print "append:", stage
        T.start()
        raw_hdf5.row_append_meas(E_l=meas[1], wtwlkr=meas[0], \
                                 proc=mea_f.proc_rank)
        tm_append += T.stop()
      #print "next"
      T.start()
      section = mea_f.seek_any_section()
      tm_seek += T.stop()
    mea_f.close()

    # Add small metadata to keep track the last successful conversion
    raw_hdf5.raw_group.attrs['convert_meas_to_hdf5::last_file'] = fname
    raw_hdf5.raw_group.attrs['convert_meas_to_hdf5::last_rank'] = mea_f.proc_rank
    rslt.flush() # maintain coherency at the end of each file

    if optTime or optDebug > 100:
      if opts.get('optimize', False):
        print "meas_count   = ", tm_meas_count(), tm_meas_count.N
      print "read_section = ", tm_readsect(), tm_readsect.N
      print "find_row     = ", tm_find(), tm_find.N
      print "append_meas  = ", tm_append(), tm_append.N
      print "seek_section = ", tm_seek(), tm_seek.N

  rslt.flush()
  return rslt


def convert_meas_to_hdf5_v2(output, H0=0, files=None, **opts):
  """A master routine to convert measurement across many pwaf-*.meas files to
  the new hdf5 format.
  Version 2, with the possibility of data reduction on-the-fly.
  Data layout is always optimized and contiguous.
  The conversion will cost a little bit more because of file opening/closing.

  The `output' parameter could be a string (interpreted as the filename of the
  target HDF5 database), or an open meas_hdf5 object.

  Additional options:
  . betablk: imaginary time length of a block
  . deltau: imaginary time step
  . debug: sets the debug level (default: 0)
  . debug_out: python file-like object to store/display the debug information
    (default: sys.stderr)
  . time: reports the timing at the end of the conversion (default: 0 == no)
  . keep_El_imag: preserves the imaginary part of E_l (default: True)
  . keep_phasefac: preserves the complex phase of wtwlkr (default: False)
  . dataset_group: the HDF5 group to contain the measurement data.
    If not specified, the default dataset group specified in the HDF5
    measurement file is used.
  . value_processor: an optional routine which takes a structured numpy array
    and cooks the values *in place* before archiving these values in the HDF5
    database. Useful for preprocessing of the values, like scaling, adding
    a constant (H0 does that), or whatnot.
    Calling convention:
      value_processor(<numpy_array>, <metadata_dict>)
    The array contains the following fields: 'E_l', 'wtwlkr', 'proc'.
    Don't mess with 'proc' unless you know what you are doing!

  Hints:
  . nwlkmax: maximum number of walkers in the measurement
  . nwlkavg: average number of walkers in the measurement
  . default_raw_chunks: a 2-tuple, i.e., '(row_chunk_size, pop_chunk_size)'
    for creating new HDF5 dataset (only used if necessary).
    measurement file is used.
    This is useful for defining a more optimal chunking parameter.

  Advanced options:
  . meas_text_class
  . meas_hdf5_class
    Used to define custom classes for creating measurement reader
    (meas_text_class) and archiver (meas_hdf5_class).
    !!! KNOW WHAT YOU ARE DOING IF YOU SPECIFY THIS !!!
    They must be API-compatible to the meas_text and meas_hdf5 defined
    in this module.
  """
  global hdf5_conv_last  # for debugging
  from pyqmc.stats.avg import avg
  from pprint import pformat

  Files = files
  files = meas_glob(Files)
  numcpu = len(files)
  if numcpu <= 0:
    raise RuntimeError, \
      "No files were found: glob = %s, pwd = %s" \
      % (str(Files), os.getcwd())

  T = timer()
  optDebug = opts.get("debug", 0)
  optTime = opts.get("time", 0)
  read_blocksize = opts.get("read_blocksize", 5)  # maximum blocks read at once
  nwlkmax = opts.get("nwlkmax", None)
  nwlkavg = opts.get("nwlkavg", None)
  debug_out = opts.get("debug_out", sys.stderr)
  keep_El_imag = opts.get("keep_El_imag", True)
  keep_phasefac = opts.get("keep_phasefac", False)
  valpx = opts.get("value_processor", None)
  meas_hdf5_class = opts.get("meas_hdf5_class", meas_hdf5)
  meas_text_class = opts.get("meas_text_class", meas_text)

  # Reads off optional parameters---these defaults better be given explicitly.
  betablk = opts.get('betablk', 0.5)
  deltau = opts.get('deltau', 0.01)

  if optDebug:
    def dbg(level, s):
      if optDebug >= level:
        debug_out.write(s)
        debug_out.flush()
  else:
    dbg = lambda level, s: None

  mea_f = meas_text_class()

  dbg(1, "Start convert_meas_to_hdf5_v2\n")
  dbg(9, "Optional parameters:\n" + pformat(opts) + "\n")

  # If nwlkmax is not provided, we will have to estimate the nwlkmax:
  # This is not a safe or exact measurement, so nwlkmax better be provided.
  nwlk_proc = []
  if nwlkmax == None:
    debug_out.write("Warning: nwlkmax not given, we will estimate it.\n")
    mea_f.open(files[0])
    nwlkmax_proc = 0
    dbg(10, "Sampling file %s...\n" % (files[0]))
    try:
      for nrecs in xrange(50):
        (stage, block, index, ndata) = mea_f.seek_any_section()
        nwlk_proc.append(ndata)
    except:
      nrecs -= 1
      pass
    nwlkmax_proc = max(nwlk_proc)
    dbg(10, "Sampled %d records, got nwlkmax_proc = %d\n" % (nrecs+1, nwlkmax_proc))
    nwlkmax = (nwlkmax_proc + 1) * numcpu * 2
    dbg(1, "Estimating nwlkmax = %d\n" % (nwlkmax))
    del nrecs

    if nwlkavg == None:
      nwlkavg = int(numpy.mean(nwlk_proc) * numcpu)
      dbg(1, "Estimating nwlkavg = %d\n" % (nwlkavg))
  else:
    if nwlkavg == None:
      nwlkavg = nwlkmax # FIXME

  if isinstance(output, meas_hdf5_class):
    rslt = output
  else:
    rslt = meas_hdf5_class()
    rslt.open(output, create_raw=False)

  # Must set the chunking parameter here if required:
  # Set chunk size to ~2% of the avg measurement count, or larger.
  default_raw_chunks = opts.get("default_raw_chunks", None)
  if default_raw_chunks != None:
    rslt.default_raw_chunks = default_raw_chunks
  dbg(10, "New chunk size = " + str(rslt.default_raw_chunks) + "\n")

  dataset_group = opts.get("dataset_group", rslt.default_raw_group)

  mea_buff = mea_f.make_numpy_buffer((read_blocksize,nwlkmax), \
                                     keep_El_imag=keep_El_imag, \
                                     keep_phasefac=keep_phasefac)

  if not dataset_group in rslt.dbh5:
    rslt.raw_create(path=dataset_group, \
                    dtypes={'E_l': mea_buff['E_l'].dtype,
                            'wtwlkr': mea_buff['wtwlkr'].dtype,
                           })
  raw_hdf5 = rslt.raw_open(path=dataset_group)
  hdf5_conv_last = rslt

  tm_readsect = avg()
  tm_find = avg()
  tm_append = avg()
  tm_seek = avg()

  read_size_actual = numpy.zeros((numcpu,), dtype=int)
  nwlkr = numpy.zeros((read_blocksize,), dtype=int)
  filepos = [None] * numcpu
  nread_total = 0
  nblk = 0
  while True:
    read_size_actual[:] = 0
    nwlkr[:] = 0
    sect_info = {}
    dbg(10, "Reading measurement records #%d\n" % (nblk))

    for icpu in xrange(len(files)):
      fname = files[icpu]
      mea_f.open(fname)
      dbg(100, "Opening file %s\n" % (fname))
      if filepos[icpu] != None:
        mea_f.file.seek(filepos[icpu])
        dbg(100, "Seeking to position %d\n" % (filepos[icpu]))

      for iblkread in xrange(read_blocksize):
        nwlk1 = nwlkr[iblkread]
        if icpu == 0:
          sect_info2 = mea_f.seek_any_section()
          sect_info[iblkread] = sect_info2
        else:
          sect_info2 = mea_f.seek_any_section()
          if sect_info2 != None:
            # check sect_info2 against sect_info[iblkread]!
            if sect_info.get(iblkread, None) == None:
              sect_info[iblkread] = sect_info2
            if sect_info[iblkread][0:3] != sect_info2[0:3]:
              sys.stderr.write("\n".join([
                "Wrong section header encountered in file " + mea_f.file.name,
                "Read header info: " + repr(sect_info2[0:3]),
                "Wanted: ", repr(sect_info[iblkread][0:3]),
                "icpu = " + str(icpu),
                "iblkread = " + str(iblkread),
              ]))
              raise RuntimeError, \
                "Wrong section header encountered in file " + mea_f.file.name

        if sect_info2 == None: break
        dbg(100, "  section: %s\n" % (str(sect_info2)))

        T.start()
        mea_f.read_section(keep_El_imag=keep_El_imag, \
                           keep_phasefac=keep_phasefac, \
                           H0=H0, out=mea_buff[iblkread, nwlk1:])
        tm_readsect += T.stop()

        dbg(100, "  read %d data\n" % (mea_f.ndata))

        mea_buff['proc'][iblkread, nwlk1 : nwlk1+mea_f.ndata] = mea_f.proc_rank

        dbg(101, "ndata,lastpos = %d %d" % (mea_f.ndata, mea_f.lastpos))
        nwlkr[iblkread] += mea_f.ndata
        filepos[icpu] = mea_f.lastpos
        read_size_actual[icpu] += 1

      mea_f.close()

    minblkread = numpy.min(read_size_actual)
    avgblkread = numpy.mean(read_size_actual)
    maxblkread = numpy.max(read_size_actual)

    if minblkread != maxblkread:
      dbg(0, \
          "At nblk = %d-%d %s:\n" % (nblk, nblk+maxblkread-1, str(sect_info[0])) + \
          "Warning: uneven number of blocks read among processors!\n" + \
          "Minimum = %d, maximum = %d, average = %.3f\n" % \
          (minblkread, maxblkread, avgblkread) + \
          "Only the minimum number of blocks are taken.\n"
          )

    for iblkread in xrange(minblkread):
      # add the measurement data to the HDF5 database
      nwlk1 = nwlkr[iblkread]
      (stage, iblk, imeas, ndata) = sect_info[iblkread]
      beta = (iblk - 1) * betablk + imeas * deltau
      extra_opts = {}
      # Creates a new measurement row
      T.start()
      row = raw_hdf5.append_row(iblk=iblk, imeas=imeas, beta=beta, \
                                stage=stage, select_row=True, \
                                **extra_opts)

      if valpx:
        meta = { # make this like HDF5 metadata for each raw data row.
          'beta': beta,
          'stage': stage,
          'iblk': iblk,
          'imeas': imeas,
          'deltau': deltau,
          'count': nwlk1,
          'row': row,
        }
        valpx(mea_buff[iblkread, 0:nwlk1], meta)

      raw_hdf5.row_append_meas(E_l=mea_buff[iblkread, 0:nwlk1]['E_l'], \
                               wtwlkr=mea_buff[iblkread, 0:nwlk1]['wtwlkr'], \
                               proc=mea_buff[iblkread, 0:nwlk1]['proc'])
      tm_append += T.stop()

    nblk += minblkread

    nread_total += minblkread
    if minblkread == 0: break
  # end main loop

  rslt.flush()

  dbg(1, "Total %d measurement blocks read\n" % (nread_total))

  if optTime or optDebug > 100:
    dbg(0, "read_section = %.3f %d\n" % (tm_readsect(), tm_readsect.N))
    dbg(0, "append_meas  = %.3f %d\n" % (tm_append(), tm_append.N))

  #raise RuntimeError,"blah"
  return rslt


def meas_glob(files, ext="ene"):
  '''Do automatic globbing on files: the files can be either a list of filenames,
  a subdir containing *.ene files, or a glob string.'''
  from glob import glob
  if files == None:
    files = sorted(glob("pwaf-*." + ext))
  elif isinstance(files, basestring):
    if os.path.isdir(files):
      files = sorted(glob(files + "/pwaf-*." + ext))
    else:
      files = sorted(glob(files))
  return files


# TENTATIVE PLACEMENT: this may be out of place:
def build_meas_index(files, **opts):
  '''Builds index of measurement data (pwaf-*.ene files).
  '''
  optDebug = opts.get("debug", 0)
  files = meas_glob(files)
  mea_f = meas_text()
  for fname in files:
    if (optDebug): print "Reading", fname
    mea_f.open(fname)
    idx_f = open(fname + ".index", "w")
    idx_f.write(\
      ["PWQMC measurement data: total energies  __INDEX_FILE__\n",
       "Process rank " + str(mea_f.proc_rank) + "\n",
       "", # UNDC HERE
      ])

    section = mea_f.seek_any_section()
    while section:
      (stage, iblk, imeas, ndata) = section
      key = ",".join( (str(stage), str(iblk), str(imeas)) )
      if optDebug > 1 and key not in n_meas_data: print "added ", key


      section = mea_f.seek_any_section()

  return n_meas_data


def get_meas_count(files, **opts):
  '''Obtains the total number of measurement data points at each
  measurement instance.
  The input is a shell wildcard of the measurement file names.
  The result can be used in convert_meas_to_hdf5 routine to increase the
  efficiency of the HDF5 archive.
  '''
  files = meas_glob(files)
  mea_f = meas_text()
  n_meas_data = {}

  optDebug = opts.get("debug", 0)

  for fname in files:
    if (optDebug): print "Reading", fname
    mea_f.open(fname)

    section = mea_f.seek_any_section()
    while section:
      (stage, iblk, imeas, ndata) = section
      key = ",".join( (str(stage), str(iblk), str(imeas)) )
      if optDebug > 1 and key not in n_meas_data: print "added ", key
      n_meas_data[key] = n_meas_data.get(key, 0) + ndata
      section = mea_f.seek_any_section()

  return n_meas_data


def print_meas_stat(meas_count):
  '''Prints a summary about measurement statistics.
  '''
  meas = numpy.array(meas_count.values())
  out = sys.stdout
  out.writelines([
    "Measurement statictics:\n",
    ". record count = %d\n" % (meas.shape[0],),
    ". nmeas avg = %d\n" % (meas.mean(),),
    ". nmeas min = %d\n" % (meas.min(),),
    ". nmeas max = %d\n" % (meas.max(),),
    ". nmeas grand total = %d\n" % (meas.sum(),),
    ". estimated HDF5 file size = %d\n" % (meas.sum() * (16+8+4),),
  ])


def reduce_measurements(wtwlkr_slices, El_slices, reblk_size):
  '''Reduces the measurement snapshots returned by load_measurements routine
  into averages, stddevs, error estimates, useful for further statistical
  analysis.'''

  from pyqmc.stats.check_reblocking import check_reblocking

  #blocks = sorted(wtwlkr_slices.keys())
  redux = {}
  reblk_sizes = [ reblk_size ]
  for blk in wtwlkr_slices:
    redux[blk] = {}
    for idx in wtwlkr_slices[blk]:
      redux[blk][idx] = check_reblocking(wtwlkr_slices[blk][idx], \
                                         El_slices[blk][idx], \
                                         reblk_sizes)[0]
  return redux


def cap_Elocal(El, El_caps, verbose = False):
  """
  Local energies should by now have been capped properly by the PWQMC-77 code.
  Should you need this functionality on a measurement dataset, please
  consider using more versatile raw_meas_data.cap_Elocal_op batch operator.
  """
  (El_min, El_max) = El_caps
  ncap_lo = 0
  ncap_hi = 0
  if verbose > 1: print "applying E_l caps from %g to %g" % (El_min, El_max)
  ndata = len(El)
  for (i, E) in zip(xrange(ndata), El):
    if E < El_min:
      ncap_lo += 1
      El[i] = El_min
    elif E > El_max:
      ncap_hi += 1
      El[i] = El_max
  if verbose and (ncap_hi > 0 or ncap_lo > 0):
    print "Total %d elements capped: %d too high, %d too low" % \
      (ncap_hi + ncap_lo, ncap_hi, ncap_lo)
  return El

