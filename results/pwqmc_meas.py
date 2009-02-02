#!/usr/bin/python
# $Id: pwqmc_meas.py,v 1.3 2009-02-02 20:02:37 wirawan Exp $
#
# pwaf_meas.py
# Tools related to PWQMC-77 measurement dump (pwaf-*.meas) files
#
# Wirawan Purwanto
# Created: 20081022
#

# Python standard modules
import sys
import time
import weakref

#import tables # HDF5-based tables--nice, but not quite...
import h5py # HDF5 Python API
import numpy

class timer:
  '''A small timer class.'''
  def start(self):
    self.tm1 = time.clock()
  def stop(self, obj):
    self.tm2 = time.clock()
    obj += (self.tm2 - self.tm1)
  def length(self):
    return self.tm2 - self.tm1


class raw_meas_rec(object):
  """Raw measurement record from one of the "load" functions.
  This class holds the basic record, which contains the array of El and
  wtwlkr data."""
  def __init__(self, El, wtwlkr):
    self.El = El
    self.wtwlkr = wtwlkr
    self.ndata = len(El)


class raw_meas_view_blk(dict):
  """View of raw measurement records based on iblk/imeas hierarchy."""
  #def __init__(self):
  def add(self, iblk, imeas, rec):
    '''Adds a measurement record to the view.'''
    self.setdefault(iblk, {})
    self[iblk][imeas] = rec
    return rec


class raw_meas_view_flat(list):
  """View of raw measurement records in a flat 1-D array."""
  def add(self, rec):
    self += [ rec ]
    return rec


# Predefined batch operations:
class batch_op(object):
  nop_global = lambda x: None  # a no-op for global level
  gather_global = lambda x: x  # a pass-through filter to get the intermediate list

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


class raw_meas_data(object):
  """Object for holding raw measurement records returned by one of
  the "load" functions.
  There are two ways of accessing the data: one is using iblk/imeas view,
  and another using flat view (no blocking)."""
  def __init__(self):
    self.blk = raw_meas_view_blk()
    self.flat = raw_meas_view_flat()
  def add(self, iblk, imeas, El, wtwlkr):
    rec = raw_meas_rec(El, wtwlkr)
    self.blk.add(iblk, imeas, rec)
    self.flat.add(rec)
    return rec
  def __getitem__(self, idx):
    '''Accesses a given measurement record. For now, it supports the
    following ways of indexing:

    object[recno]  -- accesses the "recno"-th record in the flat view.
    object[iblk,:] -- accesses the "iblk"-th block in the blk view.
    object[iblk,imeas] -- accesses the [iblk][imeas] record in the blk view.

    No slicing is supported right now.
    '''
    if type(idx) == int:
      return self.flat[idx]
    elif type(idx) == tuple:
      (iblk,imeas) = idx
      if type(imeas) == slice:
        return self.blk[iblk]
      elif type(imeas) == int:
        return self.blk[iblk][imeas]
    raise IndexError, "Invalid way of referring to a measurement record."
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
  def ndata(self):
    """Returns the number of records in the flat view."""
    return len(self.flat)
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
    Depending on the "result" setting, it can yield one of the following:
    - "array": a structured numpy array
    - "wavg":  a tuple containing weighted average of the "fld" and the sum of
      all wtwlkr"""
    if result == "array":
      global_op = \
          lambda arr: numpy.array(arr, dtype=[(fld, "=f8"), ("wtwlkr", "=f8")])
    elif result == "wavg":
      global_op = \
          lambda arr: numpy.average([ a[0] for a in arr ],
                                    weights=[ a[1] for a in arr ],
                                    returned=True)
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
    The rserve_columns option makes space for at least
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
    self.check_alive()
    E_l_ = numpy.array(E_l, dtype='=c16')

    if getattr(wtwlkr, "__iter__", False):
      wtwlkr_ = numpy.array(wtwlkr, dtype='=f8')
      if E_l_.shape != wtwlkr_.shape:
        raise TypeError, "wtwlkr array has different dimensionality from E_l"
    else:
      # Assume scalar; use broadcasting
      wtwlkr_ = numpy.empty(E_l_.shape, dtype='=f8')
      wtwlkr_[:] = wtwlkr

    if getattr(proc, "__iter__", False):
      proc_ = numpy.array(proc, dtype='=i4')
      if E_l_.shape != proc_.shape:
        raise TypeError, "proc array has different dimensionality from E_l"
    else:
      # Assume scalar
      proc_ = numpy.empty(E_l_.shape, dtype='=i4')
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
    structure in memory is also very large!'''

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
    without regard to the block/index selection.'''
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

  def load_meas_data(self, stage):
    '''Loads ALL measurement records for a given stage
    without regard to the block/index selection.
    Returns a raw_meas_data object for convenience.'''
    rslt = raw_meas_data()
    for row in xrange(0, self.row_count()):
      m = dict(zip(self.parent().raw_meta_dtype.names, self.meta[row]))
      #print m
      iblk = m['iblk']
      imeas = m['imeas']
      ndata = m['count']
      if ndata == 0: continue # don't add it in
      if m['stage'] != stage: continue
      rec = rslt.add(iblk, imeas, \
                     El=numpy.real(self.E_l[row, 0:ndata]), \
                     wtwlkr=self.wtwlkr[row, 0:ndata])
      # Also store the metadata here:
      rec.meta = m

    return rslt

class meas_hdf5(object):
  # These static fields are put here instead of in the raw_meas_hdf5
  # class, thus permitting customization of the parameters as necessary:
  raw_meta_dtype = numpy.dtype([
    ('iblk', '=i4'),
    ('imeas', '=i4'),
    ('beta', '=f8'),
    ('count', '=i4'),
    ('stage', '=i2'),
  ])
  # The default number of data is -1 (<0) as a special marker for the first
  # blank row.
  # stages: 0 = eqlb, 1 = growth, 2 = measurement
  meta_blank0 = (0, 0, 0.0, -1, 2)
  default_raw_group = 'job0/raw'
  # default chunking: 4 in time direction, 32 in wlkr direction
  default_raw_chunks = [4, 32]

  def __init__(self, fname=None, mode="a"):
    if fname: self.open(fname,mode)
  def __del__(self):
    #if hasattr(self, "dbh5")
    # We must delete circular reference to itself via raw ptr:
    #if hasattr(self, "raw"):
    #  self.raw.parent = None
    #  self.raw = None
    self.close()

  def open(self, fname, mode="a"):
    '''Opens or creates a HDF5 database containing the PWQMC measurement dump.
    This is a more efficient way of storing the pwaf-*. raw data.'''
    self.close()
    self.dbh5 = h5py.File(fname, mode)
    if "FormatType" not in self.dbh5.attrs.listnames():
      self.create_new()
    self.raw_open(self.default_raw_group)

  def close(self):
    '''Closes a previously opened HDF5 measurement database.'''
    if hasattr(self, "dbh5"):
      print "closing hdf5 database " + self.dbh5.name
      self.dbh5.close()
      delattr(self, "dbh5")
    # We must also close circular reference to itself via "raw" reference:
    if hasattr(self, "raw"):
      self.raw.parent = None
      delattr(self, "raw")

  def flush(self):
    self.dbh5.flush()

  def create_new(self):
    '''Creates a new HDF5 structure for a measurement result.'''
    db = self.dbh5
    db.attrs["Creator"] = "pyqmc.results.pwqmc_meas"
    db.attrs["CreateTime"] = \
        time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
    db.attrs["FormatVersion"] = 0
    db.attrs["FormatType"] = "pwqmc_meas"
    db.require_group("job0")
    self.raw_create(self.default_raw_group)
    db.flush()

  def raw_create(self, path):
    '''Creates a new raw measurement section at absolute group location `path'.'''
    g = self.dbh5.require_group(path)
    '''
    Notes for arrays E_l, wtwlkr, and proc:
    * these all must have identical dimensions
    * index #1 `meas_idx': global measurement index (0, 1, 2, 3, ...)
    * index #2 `wlkr_idx': walker index (0, 1, 2, 3, ..., meta[meas_idx].count)

    The meta array runs with meas_idx only.
    '''

    Chunk = tuple(self.default_raw_chunks)
    g.create_dataset("E_l", shape=(1, 1), dtype="=c16", \
        maxshape=(None,None), chunks=Chunk)
    g.create_dataset("wtwlkr", shape=(1, 1), dtype="=f8", \
        maxshape=(None,None), chunks=Chunk)
    g.create_dataset("proc", shape=(1, 1), dtype="=i4", \
        maxshape=(None,None), chunks=Chunk)
    # Measurement markers (blk/imeas/ compound):
    meta_blank = self.raw_meta_blank()[0]
    g.create_dataset("meta", shape=(1,), dtype=self.raw_meta_dtype, \
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
    self.raw = raw  # last opened raw group (WARNING: circular reference)
    return raw

  


# end class meas_hdf5



class meas_text(object):
  '''meas_text is a class for reading the measurement dump produced by
  PWQMC-77 code (i.e. a single pwaf-*.meas or pwaf-*.ene file).
  The raw dump is in text format.'''
  def __init__(self, fname = None):
    if fname: self.open(fname)
  def __del__(self):
    self.close()

  # Opening and closing input files:
  def open(self, fname):
    '''Opens a given measurement file for reading and analysis.'''
    self.close()
    self.file = open(fname, "r")
    #self.fname = fname  -- use self.file.name instead
    # Requires a proper header:
    F = self.file
    L = F.readline().strip()
    if not L.startswith("PWQMC measurement data: total energies"):
      raise ValueError, "File `" + fname + "': Invalid file type"

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

  def close(self):
    '''Closes a previously opened measurement file.'''
    if hasattr(self, "file"):
      self.file.close()
      delattr(self, "file")

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
          return True
      line = F.readline()
    return False

  def seek_any_section(self):
    '''Seeks forward the any next named section.
    Returns a tuple of 4 elements: (stage, block, index, ndata).
    Stage is an integer of value 0 (equilibration), 1 (growth), or 2
    (measurement). Block is the block number (0, 1, 2, ...).
    Index indicates the measurement index. (??? <<< FIXME DOC HERE)
    Ndata is the number of data points in this measurement block (per
    MPI process).
    On the earliest (obsolete) dataset format, Ndata field is nonexistent,
    and it is set to None.
    Returns None if there is no more measurement section detected.
    '''
    F = self.file
    line = F.readline()
    #for line in self.file:
    while line != "":
      # This kind of seek should be very fast:
      L = line.lstrip()
      if L.startswith('MEAS ') or L.startswith('EQLB ') or L.startswith('GRTH '):
        flds = line.split()
        if L.startswith('MEAS '):
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
        return (stage, block, index, ndata)
      line = F.readline()
    return None

  def read_section(self, section = "MEA", block = 1, index = 1, \
    fname = None, keep_El_imag = False, H0 = 0):
    '''Reads energy measurement, one section at a time from a pwaf-*.meas file.
    Returns a 2-component tuple containing (1) wtwlkr and (2) local energy
    data.

    After the end of this subroutine, the file pointer is positioned to read
    the next section right after it.

    If section is set to None, then it begins reading the measurement section
    immediately.
    '''
    if fname != None: self.open(fname)
    if section != None and not self.seek_section(section, block, index):
      raise ValueError, \
        "Cannot seek section " + section + " block " \
        + str(block) + " index " + str(index)

    wtwlkr = []
    E_l = []
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
      wtwlkr.append(float(flds[1]))
      if (keep_El_imag):
        E_l.append(complex(float(flds[2]), float(flds[3])) + H0)
      else:
        E_l.append(float(flds[2]) + H0)
      lastpos = F.tell()
      line = F.readline()

    rslt = (wtwlkr, E_l)
    return rslt


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
  memory is also very large!'''
  from glob import glob

  if files == None:
    files = sorted(glob("pwaf-*.meas"))
  elif type(files) == str:
    files = sorted(glob(files))

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


def convert_meas_to_hdf5(output, H0 = 0, files = None, **opts):
  '''A master routine to convert measurement across many pwaf-*.meas files to
  the new hdf5 format.

  Additional options:
  . betablk: imaginary time length of a block
  . deltau: imaginary time step
  . optimize: finds out the number of measurements in all files
    and use that information to optimize the storage use and layout.
    This can be a very long operation if the files are large or many.
  . nmeas_counts: provides the number of measurement points in each record.
    Usually this is not necessary as it can be provided via optimize()
    above.
  '''
  global hdf5_conv_last
  from glob import glob
  from pyqmc.stats.avg import avg

  if files == None:
    files = sorted(glob("pwaf-*.meas"))
  elif type(files) == str:
    files = sorted(glob(files))

  T = timer()
  tm_meas_count = avg()

  optDebug = opts.get("debug", 0)
  optTime = opts.get("time", 0)
  if opts.get('optimize', False):
    T.start()
    nmeas_counts = get_meas_count(files)
    T.stop(tm_meas_count)
  else:
    nmeas_counts = opts.get("nmeas_counts", None)

  mea_f = meas_text()
  rslt = meas_hdf5()
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
      #print "read_sect"
      T.start()
      meas = mea_f.read_section(section=None, H0=H0, keep_El_imag=True)
      T.stop(tm_readsect)
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

      T.stop(tm_find)
      if len(meas[0]) > 0:
        #print "append:", stage
        T.start()
        raw_hdf5.row_append_meas(E_l=meas[1], wtwlkr=meas[0], \
                                 proc=mea_f.proc_rank)
        T.stop(tm_append)
      #print "next"
      T.start()
      section = mea_f.seek_any_section()
      T.stop(tm_seek)
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


# TENTATIVE PLACEMENT: this may be out of place:
def get_meas_count(files, **opts):
  '''Obtains the total number of measurement data points at each
  measurement instance.
  The input is a shell wildcard of the measurement file names.
  The result can be used in convert_meas_to_hdf5 routine to increase the
  efficiency of the HDF5 archive.
  '''
  from glob import glob

  if files == None:
    files = sorted(glob("pwaf-*.meas"))
  elif type(files) == str:
    files = sorted(glob(files))

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

