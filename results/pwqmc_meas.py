#!/usr/bin/python
# $Id: pwqmc_meas.py,v 1.2 2009-01-21 21:23:45 wirawan Exp $
#
# pwaf_meas.py
# Tools related to PWQMC-77 measurement dump (pwaf-*.meas) files
#
# Wirawan Purwanto
# Created: 20081022
#

# Python standard modules
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

class raw_meas_hdf5(object):
  """Creates a new accessor for a `raw measurement' data structure.
  This object is not to be created directly by use; instead, use
  meas_hdf5.raw_open() object method to create the accessor.
 
  NOTE: the accessor object's validity depends on the life of the parent
  meas_hdf5 object."""
  def __init__(self, raw_group, parent):
    """Creates a new raw measurement accessor.
    This routine is to be called by the meas_hdf5.raw_open() object method.
    """
    # The reference to parent *MUST* be a weak reference or else it would
    # cause the
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

  def append_row(self, iblk, imeas, beta, stage = 2, select_row = False):
    '''Appends a new row on the measurement record.
    Optionally, selects the new row as the current row.
    Returns the integer index of the new measurement row.'''
    self.check_alive()
    g = self.raw_group
    meta_new = numpy.array((iblk, imeas, beta, 0, stage), \
                           dtype=self.parent().raw_meta_dtype)

    do_extend = True
    shape = self.E_l.shape
    if shape[0] == 1: # may contain a dummy empty element
      meta = self.meta[0]
      if meta[3] < 0:
        do_extend = False
        self.meta[0] = meta_new

    if do_extend:
      shape = self.extend_shape((shape[0]+1, shape[1]))
      #self.E_l.id.extend(shape)
      #self.wtwlkr.id.extend(shape)
      #self.proc.id.extend(shape)
      #self.meta.id.extend(shape[0:1])
      self.meta[shape[0]-1] = meta_new

    if select_row: self.row = shape[0]-1
    return shape[0]-1

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
    E_l_ = numpy.array(E_l, dtype='=c8')

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

  def load_meas(self, blocks, indices, H0 = 0, stage = 2):
    '''A master routine to load measurement across many blocks/indices
    into memory (compatible with early "slices" format).
    Returns the snapshots of (1) wtwlkr and (2) local energy as a tuple.
  
    Be aware that if the measurement files are large, then the data
    structure in memory is also very large!'''

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
      print "closing hdf5 database ", self.dbh5.name
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
    g.create_dataset("E_l", shape=(1, 1), dtype="=c8", \
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
    '''Seeks forward the any next named section.'''
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
        return (stage, block, index)
      line = F.readline()
    return False

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
  '''
  global yyy
  from glob import glob
  from pyqmc.stats.avg import avg

  if files == None:
    files = sorted(glob("pwaf-*.meas"))
  elif type(files) == str:
    files = sorted(glob(files))

  mea_f = meas_text()
  rslt = meas_hdf5()
  rslt.open(output)
  raw_hdf5 = rslt.raw
  cache = raw_hdf5.build_row_cache('stage', 'iblk', 'imeas')
  yyy = rslt
  betablk = opts.get('betablk', 0.5)
  deltau = opts.get('deltau', 0.01)

  tm_readsect = avg()
  tm_find = avg()
  tm_append = avg()
  tm_seek = avg()
  T = timer()
  for fname in files:
    #if f > "pwaf-00001.meas": break
    print "Reading", fname
    mea_f.open(fname)
    print "file ", fname, "rank = ", mea_f.proc_rank, "start @ ", mea_f.start_date, mea_f.start_time
    #print "<" + line + ">"

    section = mea_f.seek_any_section()
    while section:
      (stage, iblk, imeas) = section
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
      '''
      row = raw_hdf5.find_row(iblk=iblk, imeas=imeas, stage=stage, \
                              _create_new = {'beta': beta}, \
                              _select_row = True)
      '''
      # For now we manage the cache ourselves. Maybe later we can use a
      # different means.
      key = ",".join( (str(stage), str(iblk), str(imeas)) )
      if key not in cache:
        row = raw_hdf5.append_row(iblk=iblk, imeas=imeas, beta=beta, \
                                  stage=stage, select_row=True)
        cache[key] = row
      else:
        raw_hdf5.select_row(cache[key])

      T.stop(tm_find)
      if len(meas[0]) > 0:
        #print "append:", stage
        T.start()
        raw_hdf5.row_append_meas(E_l=meas[1], wtwlkr=meas[0], proc=mea_f.proc_rank)
        T.stop(tm_append)
      #print "next"
      T.start()
      section = mea_f.seek_any_section()
      T.stop(tm_seek)
    mea_f.close()

    print "read_section = ", tm_readsect(), tm_readsect.N
    print "find_row     = ", tm_find(), tm_find.N
    print "append_meas  = ", tm_append(), tm_append.N
    print "seek_section = ", tm_seek(), tm_seek.N
  return rslt


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
  if verbose: print "applying E_l caps from %g to %g" % (El_min, El_max)
  ndata = len(El)
  for (i, El) in zip(xrange(ndata), El):
    if El < El_min:
      ncap_lo += 1
      El_w[i] = El_min * wtwlkr[i]
    elif El > El_max:
      ncap_hi += 1
      El_w[i] = El_max * wtwlkr[i]
  if verbose and (ncap_hi > 0 or ncap_lo > 0):
    print "Total %d elements capped: %d too high, %d too low" % \
      (ncap_hi + ncap_lo, ncap_hi, ncap_lo)


