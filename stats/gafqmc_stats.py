#
# pyqmc.stats.gafqmc_info
# Tools to parse GAFQMC stats file
#
# Wirawan Purwanto
# Created: <200911 in Cr2_BFD_analyze1.py
# Imported: 20130920 into pyqmc.
#

"""
pyqmc.stats.gafqmc_info
Tools to parse and process GAFQMC stats file.
"""

import numpy
import re
import sys

from pyqmc import PyqmcParseError
from wpylib.file.file_utils import open_input_file

class gafqmc_stat(object):
  '''GAFQMC statistics file.'''
  # Timing fields which are averaged over time:
  # These are all the possible timing fields in the stats file.
  timing_field_map = {
    '<step>': 'step',
    '<overlap>': 'ovlp',
    '<Calc_FB>': 'FB',     # renamed (CalcFB -> FB)
    '<Calc_HSop>': 'HSop', # not existing in fort.17
    '<ExpH>': 'ExpH',
    '<Elocal>': 'Elocal',
    '<gasdev_auxfld>': 'gasdev',
    # MPI version:
    '<pctl>': 'pctl',
    '<balance>': 'pc_balance',
    '<comm:slave2master>': 'pc_s2m',
    '<comm:master2slave>': 'pc_m2s',
    '<comm:slave2slave>': 'pc_s2s',
    '<comm:barrier>': 'pc_barrier1',
    # Newer MPI version >= 20131202 #2486f9fe have additional fields
    '<comm:barrier2>': 'pc_barrier2',
    '<modGS>': 'modGS',
    '<AccBlk>': 'AccBlk',
    '<write_chkpt>': 'chkpt',
  }

  #rx_table_header = re.compile(r'^\s*(Equilibration|Growth|Measurement) phase\.\.\.')
  rx_table_header = re.compile(r'^\s*(Equilibration|Growth|Measurement) phase\.\.\.|^\s*t_marker: (Equilibration|Growth|Measurement)_phase at ')

  debug = 10
  def dbg(self, level, msg):
    if self.debug >= level:
      sys.stderr.write(msg + "\n")

  def read(self, filename, save_text=False):
    self.save_text = save_text
    if save_text:
      self.text = {}
    else:
      if hasattr(self, "text"):
        delattr(self, "text")
    self.stats_blk = {}  # only exists if decumulate is "True" (default)
    self.stats_acc = {}
    F = open_input_file(filename)

    self.Parse_preamble_(F)

    while True:
      if self.parse_last_text_ == None:
        break
      elif self.parse_phase_ == "POSTAMBLE":
        self.Parse_postamble_(F)
        break
      elif self.parse_phase_ in ("Equilibration", "Growth", "Measurement"):
        r = self.Parse_table_(F)
      else:
        raise PyqmcParseError, "BUG: Unknown parse_phase = %s" % self.parse_phase_

    F.close()

  def Parse_detect_table_header_(self, text):
    m = self.rx_table_header.search(text)
    if m:
      self.dbg(2, "Matched table header text: %s" % text.rstrip())
      self.dbg(10, "m.groups = %s" % (m.groups(),))
      for m1 in m.groups():
        if m1 != None:
          self.parse_next_phase_ = m1
          return True
      raise PyqmcParseError, "BUG: cannot find next phase?"
    else:
      return False

  def Parse_preamble_(self, F):
    save_text = self.save_text
    self.parse_phase_ = "PREAMBLE"
    if save_text:
      text = []
      self.text[self.parse_phase_] = text

    for L in F:
      L = L.rstrip()
      if self.Parse_detect_table_header_(L):
        self.parse_last_text_ = L
        self.parse_phase_ = self.parse_next_phase_
        return self.parse_phase_
      if save_text:
        text.append(L)
      # TODO: do file preamble parsing here

    self.parse_last_text_ = None

  def Parse_postamble_(self, F):
    save_text = self.save_text
    if self.save_text:
      text = []
      self.text[self.parse_phase_] = text
      if self.parse_last_text_ != None:
        text.append(self.parse_last_text_)
      for L in F:
        L = L.rstrip()
        if save_text:
          text.append(L)
        # TODO: do file postamble parsing here

  def Parse_table_(self, F):
    """Parses a stats main table after the preamble or another table.
    """
    save_text = self.save_text
    phase = self.parse_phase_
    self.dbg(1, "Parsing table: %s" % phase)
    rx_tbl_hdr = re.compile(r'Block\s+Time_elapsed')
    tbl_dtype = None
    if save_text:
      text = []
      self.text[phase] = text
      if self.parse_last_text_ != None:
        text.append(self.parse_last_text_)

    self.parse_last_text_ = None
    for L in F:
      # TODO: do table preamble parsing here
      Ls = L.strip()
      self.dbg(100, "tbl: "+Ls)
      if Ls.startswith("timing: Statistical snapshots"):
        pass
      elif Ls.startswith("t_marker: "+phase+" phase at") \
           or Ls.startswith(phase+" phase..."):
        pass
      elif Ls.startswith("timing: "+phase+" phase ="):
        pass
      elif self.Parse_detect_table_header_(L):
        # set up for next parsing phase
        self.parse_last_text_ = L.rstrip()
        self.parse_phase_ = self.parse_next_phase_
        return -1
      elif rx_tbl_hdr.match(Ls):
        tbl_dtype = self.make_table_fields(Ls)
      elif Ls[0].isdigit():
        # Take this as the first text of the table
        self.parse_last_text_ = L.rstrip()
        if tbl_dtype == None:
          raise PyqmcParseError, \
            "Cannot determine table structure for phase: %s" % (phase,)
        self.Parse_table_content_(F, tbl_dtype)
        self.parse_phase_ = None
        # Next table, maybe?
        if self.parse_last_text_ != None:
          if self.Parse_detect_table_header_(self.parse_last_text_):
            # set up for next parsing phase
            self.parse_phase_ = self.parse_next_phase_
          else:
            # Otherwise, assume a postamble afterward
            self.parse_phase_ = "POSTAMBLE"
        return +1
      else:
        # Other things
        self.parse_last_text_ = L.rstrip()
        self.parse_phase_ = "POSTAMBLE"
        break
        #raise PyqmcParseError, \
        #  "Invalid table text in phase '%s': %s" % (self.parse_phase_,Ls)

      if save_text:
        self.dbg(100, "add_text: %s %s" % (phase, L.rstrip()))
        text.append(L.rstrip())

    # If you reach this point, something might have been parsed
    # but not the table's content
    return 0

  def Parse_table_content_(self, F, dtype, decumulate=True):
    """Parses and loads a single stats table into an array.
    The output records are still cumulative in nature unless
    decumulate is True.

    NOTE: To save space, we do NOT save the text of this section.
    """
    import itertools
    phase = self.parse_phase_
    rslt0 = []
    # Remove all NaNs and return them as zero
    NaN_to_zero = lambda x : x if not numpy.isnan(x) else 0
    if getattr(self, "parse_last_text_", None) != None:
      F0 = [ self.parse_last_text_ ]
    else:
      F0 = []
    self.parse_last_text_ = None
    for ll in itertools.chain(F0, F):
      L = ll.strip()
      # Right now simply use this criterion as a marker of a record, i.e.
      # the first field being a number: it's not foolproof, man:
      if L[0].isdigit():
        rslt0.append(tuple(L.split()))
      else:
        self.parse_last_text_ = ll.rstrip()
        break

    if rslt0:
      self.stats_acc[phase] = numpy.array(rslt0, dtype=dtype)
      if decumulate:
        self.stats_blk[phase] = stat_decumulate(self.stats_acc[phase])
      self.dbg(2, "Parsed table for phase: %s, %d rows" % (phase, len(self.stats_acc[phase])))
    else:
      self.dbg(2, "Warning: empty table for phase: %s" % phase)
      return None


  def make_table_fields(self, header_text):
    """Scans the table header text and constructs the right datatype for
    this table.
    """
    flds = header_text.split()
    if len(flds) > 2 and flds[0] == "Block" and flds[1] == "Time_elapsed":
      # Fixed portion
      stat_dtype0 = [
        ('idx', int),
        ('T_elapsed', float),
      ]
    else:
      raise PyqmcParseError, \
        "Invalid gafqmc_stats table header prefix: %s %s" % (flds[0], flds[1])
    # Scan the variable timing fields
    self.timing_fields = []
    for f in flds[2:]:
      try:
        K = self.timing_field_map[f]
      except KeyError:
        raise PyqmcParseError, \
          "Invalid gafqmc_stats table field name: %s" % (f,)
      self.timing_fields.append(K)

    if self.timing_fields[-1] == "pc_barrier1":
      # Older MPI GAFQMC have this at the end, but actually there are 2 barrier fields:
      self.timing_fields.append("pc_barrier2")

    self.dbg(2, "  Fields=%d %s" % (len(self.timing_fields), " ".join(self.timing_fields)))

    for f in self.timing_fields:
      stat_dtype0 += [ ('t_' + f, float), ('n_' + f, int) ]

    self.stat_dtype = numpy.dtype(stat_dtype0)
    return self.stat_dtype


def stat_decumulate(rec):
  """Converts the stat record from the cumulative (original) to
  per-block quantities.
  Returns a copy object containing the decumulated quantities."""
  rslt = rec.copy()
  #print rslt
  t_elapsed = rslt['T_elapsed']
  # Unfortunately numpy does not have a converse function of "accumulate".
  # Work on them in batch as much as possible:
  x_timing_fields = rec.dtype.names[2:]
  Nfields = len(x_timing_fields) // 2
  tval = numpy.empty((Nfields, len(rec)), dtype=float)
  nval = numpy.empty((Nfields, len(rec)), dtype=int)
  for (i,col_t, col_n) in zip(xrange(Nfields), x_timing_fields[::2], x_timing_fields[1::2]):
    tval[i] = rec[col_t]
    nval[i] = rec[col_n]
  #print tval; print nval
  # make tval a sum instead of the accumulated average
  tval *= nval
  #print tval
  # "de"-cumulate the sums
  for i in xrange(len(rec)-1, 0, -1):
    tval[:,i] -= tval[:,i-1]
    nval[:,i] -= nval[:,i-1]
    t_elapsed[i] -= t_elapsed[i-1]
  nval_nz = nval.copy()
  #print tval; print nval
  # Avoid division by zero:
  numpy.putmask(nval_nz, nval_nz==0, 1)
  # Zero time averages where no events recorded:
  numpy.putmask(tval, nval_nz==0, 0.0)
  # reconstruct the block averages
  tval /= nval_nz
  #print nval_nz
  for (i,col_t, col_n) in zip(xrange(Nfields), x_timing_fields[::2], x_timing_fields[1::2]):
    rslt[col_t] = tval[i]
    rslt[col_n] = nval[i]
  return rslt

