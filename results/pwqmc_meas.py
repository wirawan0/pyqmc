#!/usr/bin/python
# $Id: pwqmc_meas.py,v 1.1 2009-01-09 22:06:14 wirawan Exp $
#
# pwaf_meas.py
# Tools related to PWQMC-77 measurement dump (pwaf-*.meas) files
#
# Wirawan Purwanto
# Created: 20081022
#

import numpy

class pwqmc_meas_file(object):
  '''pwqmc_meas_file is a class for reading the measurement dump produced by
  PWQMC-77 code (i.e. a single pwaf-*.meas file).'''
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

  def read_section(self, section = "MEA", block = 1, index = 1, \
    fname = None, keep_El_imag = False, H0 = 0):
    '''Reads energy measurement, one section at a time from a pwaf-*.meas file.
    Returns a 2-component tuple containing (1) wtwlkr and (2) local energy
    data.

    After the end of this subroutine, the file pointer is positioned to read
    the next section right after it.
    '''
    if fname != None: self.open(fname)
    if not self.seek_section(section, block, index):
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
  if not (idx in dest):
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
  mea_f = pwqmc_meas_file()
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


