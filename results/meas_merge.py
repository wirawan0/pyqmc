#
# pyqmc.results.meas_merge
# Tools to merge raw measurements found in AFQMC INFO file
#
# Wirawan Purwanto
# Created: 20130226
#

"""\
pyqmc.results.meas_merge
Utilities to merge raw measurements found in AFQMC INFO file.
"""

import numpy
import sys

from wpylib.sugar import sorted

class merged_measurement(object):
  """Object to hold merged measurement results.
  """
  rdiff = 1e-9
  # Similar to my get-energy.gawk --merge
  # Custom meas_dtype like that from gafqmc_info, with extra fields
  meas_dtype = numpy.dtype([('beta',float), ('overlap',float),
                            ('Etotal',float), ('Eproj',float),
                            ('ndata',int)])
  default_fix_format = "%.12f"

  def merge(self, info_recs):
    """Merges raw measurement records found in 'INFO' structures
    (gafqmc_info, pwqmc_info).
    """
    self.beta_merged_set = set()
    for info in info_recs:
      self.collect_beta(info)
    self.consolidate_beta(rdiff=self.rdiff)
    self.merged = numpy.zeros((len(self.betas),), dtype=self.meas_dtype)
    self.merged['beta'] = self.betas

    bmap = self.beta_mapping
    merge = self.merged
    for info in info_recs:
      for R in info:
        (beta, overlap, Etotal, Eproj) = (R['beta'], R['overlap'], R['Etotal'], R['Eproj'])
        i = bmap[beta]
        D = merge[i]
        D['overlap'] += overlap
        D['Etotal'] += Etotal * overlap
        D['Eproj'] += Eproj * overlap
        D['ndata'] += 1
    for D in merge:
      overlap = D['overlap']
      if overlap:
        D['Etotal'] /= overlap
        D['Eproj'] /= overlap

    return merge

  def collect_beta(self, raw_rec):
    """Collect the beta values (unedited) to the collection bin
    (self.beta_merged_raw)."""
    beta = raw_rec['beta']
    bmerge = self.beta_merged_set
    bmerge.update(set(beta))

  def consolidate_beta(self, rdiff):
    """Consolidates the beta values so that nearby values are fused."""
    from wpylib.db.indexing_float import generate_float_indices
    rec = generate_float_indices([ b for b in self.beta_merged_set ], rdiff_threshold=rdiff)
    self.betas = rec.vals
    self.beta_mapping = rec.index_mapping

  def print_formatted(self, fmt, out=None):
    if out == None: out = sys.stdout
    for D in self.merged:
      out.write(fmt % (D))
    try:
      Flush = out.flush
    except:
      pass
    else:
      Flush()


def autodetect_fix_format(self, values, fmt):
  """Tries to find the format that would display all the values in a nice
  columnar fashion.

  Warning: this routine could be very slow or use lots of memory for a very
  large dataset!
  """
  svals = [ (fmt % v).split(".",1) for v in values ]







