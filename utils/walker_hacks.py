#
# pyqmc.utils.walker_hacks module
#
# Wirawan Purwanto
# Created: 20130625
#

"""
pyqmc.utils.walker_hacks

Hacks for walker files and data.
"""

import os
import os.path
import sys

import numpy.random

from glob import glob
from os.path import abspath, basename, dirname
from os.path import join as joinpath

import wpylib.shell_tools as sh

import wpylib.iofmt.fortbin
from wpylib.iofmt.text_output import text_output


ftypes = wpylib.iofmt.fortbin.fortran_types()
# tweak ftypes to adjust your needs (e.g. if ILP64 system is in effect).

class emergency_walkers_fixup(object):
  """A tool to quickly fix a collection of QMC walkers for reuse in a subsequent
  QMC run.

  Motivation: In very large computers sometimes walker files can get corrupted.
  In this case, to prevent the QMC job from failing, we simply 'sample' the
  walker population to regenerate the missing/broken walker file.

  This tool provides a general framework to perform this emergency saving
  in the case of broken files.
  """

  # walkers_fixup_dir is a relative subdirectory to contain fixed up walker file set
  # walkers_fixup_note_bad_dir is a relative subdirectory to contain ONLY the replaced
  # bad walker filenames (for human inspection)
  walkers_fixup_dir = "fixup"
  walkers_fixup_note_bad_dir = "fixup/bad"

  verbose = 1
  log = text_output(sys.stdout, flush=True)

  # The following must be overriden:
  walkers_glob = None

  """
  Transient variables:
  - walkers_dir
  - walkers_absdir
  - walkers[]
  - walkers_bad[]
  - walkers_ibad[] -- index of bad walkers
  - walkers_good[]
  - walkers_reserved[] -- names of walker files that are NOT to be used
    as replacement for bad walkers.
  """

  def is_bad_walker_file(self, filename):
    """Decides if a walker file is bad. Must be overriden in actual class."""
    return 0

  def collect_walker_files(self, dir="wlk"):
    orig_pwd = os.getcwd()
    self.walkers_absdir = abspath(dir)
    os.chdir(dir)
    self.walkers_dir = dir
    self.walkers = sh.sorted_glob(self.walkers_glob)
    if self.verbose:
      self.log(": Found %d walker files\n" % len(self.walkers))
    os.chdir(orig_pwd)
    try:
      del self.walkers_bad
      del self.walkers_ibad
      del self.walkers_good
    except:
      pass

  def check_walker_files(self, report=None):
    """Checks the walker files (can take awhile), and returns the number of
    bad walker files found."""
    self.walkers_bad = []
    self.walkers_ibad = []
    self.walkers_good = []
    for (iw,w) in enumerate(self.walkers):
      if self.is_bad_walker_file(w):
        self.walkers_bad.append(w)
        self.walkers_ibad.append(iw)
      else:
        self.walkers_good.append(w)
    return len(self.walkers_bad)

  @property
  def has_bad_walkers(self):
    if not hasattr(self, "walkers_bad"):
      self.check_walkers_bad(report=False)
    return len(self.walkers_bad) > 0

  def report_bad_walkers(self, level):
    if self.has_bad_walkers:
      self.log("! Found %d bad walker files...\n" % len(self.walkers_bad))
      if level > 0 and level < 10:
        self.log("".join(["! Bad files (indices): "]
                         + [" %s" % iw for iw in self.walkers_ibad ] + ["\n"]))
      else:
        self.log("".join(["! Bad files:\n"]
                         + ["! - %6d  %s\n" % (iw,w) 
                              for (iw,w) in zip(self.walkers_ibad, self.walkers_bad) ]))

  def fixup_walkers1(self):
    #if not self.has_bad_walkers:
    #  return False

    orig_pwd = os.getcwd()
    os.chdir(self.walkers_absdir)
    sh.mkdir("-p", self.walkers_fixup_dir)
    sh.mkdir("-p", self.walkers_fixup_note_bad_dir)

    # Copy over the good walkers
    for w in self.walkers_good:
      self.copy_file(w, joinpath(self.walkers_fixup_dir,w))

    num_good = len(self.walkers_good)
    max_tries = 100

    for wb in self.walkers_bad:
      tries = 0
      while tries <= max_tries:
        tries += 1
        r = numpy.random.randint(num_good)
        wg = self.walkers_good[r]
        if wg not in self.walkers_reserved:
          break
        elif tries == max_tries:
          raise RuntimeError, \
            "Fatal: cannot find good replacement for walker file %s" % wb
      if self.verbose >= 10:
        self.log("fixup: %s -> %s\n" % (wb, wg))
      self.copy_file(wg, joinpath(self.walkers_fixup_dir, wb))
      self.copy_file(wg, joinpath(self.walkers_fixup_note_bad_dir, wb))

    os.chdir(orig_pwd)

  def generate_replacement_list(self):
    """Generates a list of walker replacement. Excluding those that are in
    walkers_reserved array.
    """
    raise NotImplementedError


  def copy_file(self, src, dest):
    from wpylib.file.file_utils import relpath
    # FIXME--may want real copy command, perhaps?
    sh.provide_link(dest, relpath(dirname(abspath(dest)), abspath(src)))



class gafqmc_emergency_walkers_fixup(emergency_walkers_fixup):
  """Emergency walker population fixer for GAFQMC code."""
  """
  FIXME:
   - filenames for >= 1e6 processors
  """
  walkers_glob = 'gafqmc-' + ("[0-9]" * 5)
  walkers_reserved = [ 'gafqmc-00000' ]

  # Structure of walker file header, processes 1..Nproc-1
  struct_hdr_walker_file_slave = [
    (int,int,int,int), # lran1..4
    (float,), # version_chkpt
    (int, float, float, int), # nh, anorm, etrial, istpacc
    (float,), # times
    (int,), # counts
  ]
  def calc_min_walker_file_size(self):
    return ftypes.file_data_size(self.struct_hdr_walker_file_slave)

  @property
  def min_walker_file_size(self):
    if not hasattr(self, "_min_walker_file_size"):
      self._min_walker_file_size = self.calc_min_walker_file_size()
    return self._min_walker_file_size

  def is_bad_walker_file(self, filename):
    """Decides if a walker file is bad.
    For now we decide this merely based on file size."""
    return os.stat(joinpath(self.walkers_absdir, filename)).st_size < self.min_walker_file_size

