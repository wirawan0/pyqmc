# -*- python -*-
#
# pyqmc.utils.gafqmc_dir.py
# Standard convention of GAFQMC calculations
#
# Wirawan Purwanto
# Created: 20140115
#
# This module is part of PyQMC, a toolsuite for AFQMC projects.

"""pyqmc.utils.gafqmc_dir.py
Standard convention of GAFQMC calculations

This module is part of PyQMC, a toolsuite for AFQMC projects.

This module defines the author's standard convention
for GAFQMC's file organization.
Adherence to this convention is strongly recommended so you can use
all the tools available in the pyqmc python package.

"""

# Python standard modules
import os
import os.path
import re
import sys
import time

import numpy

import wpylib.shell_tools as sh
from pyqmc import PyqmcWarning

gafqmc_out_file_patterns = set(('INFO','INFO.gz','INFO.bz2','INFO.lzma','INFO.xz'))


def is_gafqmc_result_dir(D, files=None, dirs=None,
                         file_pattern=None, parse_file=True):
  """Tests whether the directory D, containing `files' (and softlinks)
  and directories `dirs' is a result directory for a GAFQMC-type
  calculation.
  Returns the score of the test, where higher score means more
  reliability.

  Input arguments `files' and `dirs' are not necessary (in fact, not
  recommended) unless you use this in conjunction with os.walk,
  where the files and dirs would have been gathered during the
  iteration cycle.

  Return flag: an integer or-ed
    1 = output file name exists
    2 = AND output file does exists as a regular file
    4 = AND output file is indeed a GAFQMC output file
    8 = AND input filename exists that matches the output
  """
  from os.path import join, isdir, isfile
  from wpylib.sugar import is_iterable
  from wpylib.file.file_utils import list_dir_entries
  from pyqmc.results.gafqmc_info import is_gafqmc_info

  if files == None or dirs == None:
    dirs, files = list_dir_entries(D)[:2]

  rslt = 0
  if file_pattern == None:
    file_pattern = gafqmc_out_file_patterns
  if isinstance(file_pattern, (set, tuple, list)) or is_iterable(file_pattern):
    if not isinstance(file_pattern, set):
      file_pattern = set(file_pattern)
    fset = set(files)
    fset_good = file_pattern & fset
    if len(fset_good) > 0:
      # WARNING: This will create uncertainty if there are more than one file
      # matching the pattern. BE WARNED!
      info_file = sorted(list(fset_good))[0]
      rslt |= 1
  else:
    raise NotImplementedError

  if rslt:
    # At least the filename is found:
    info_path = join(D, info_file)
    if isfile(info_path):
      rslt |= 2
      if parse_file and is_gafqmc_info(info_path):
        rslt |= 4
        # the next if's are TO BE IMPLEMENTED LATER

  return rslt


def get_latest_rundir(rsltdir):
  """Obtains the latest rundir based on mtime of the rundir symlinks.
  The input is a QMC result dir.
  """
  from warnings import warn
  from glob import glob
  from os.path import join, islink, isdir
  rsltpath = lambda x : join(rsltdir, x)
  def cmp_mtime(r1, r2):
    return cmp(r1[1], r2[1])
  candidates = []
  for p in [ rsltpath("rundir") ] + glob(rsltpath(".rundir.[0-9]*.[0-9]*")):
    if islink(p) or isdir(p):
      try:
        mtime = os.lstat(p).st_mtime
      except OSError:
        continue
      candidates.append((p, mtime))
  if len(candidates) == 0:
    return None
  candidates.sort(cmp=cmp_mtime, reverse=1)
  if len(candidates) > 1:
    sys.stderr.write("get_latest_rundir: found %d rundirs, choosing the most recent one:\n -> %s\n" \
                     % (len(candidates), "\n -  ".join([ s[0] for s in candidates ])))
    warn("get_latest_rundir: multiple rundirs detected", PyqmcWarning)
  return candidates[0][0]



