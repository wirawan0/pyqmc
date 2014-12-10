#!/usr/bin/python
#
# pyqmc.sites.sciclone.tools module
#
# Wirawan Purwanto
# Created: 20110925
#

"""
pyqmc.sites.sciclone.tools

Site-specific support for pyqmc library.
This module is part of PyQMC project.

This module provides specific tools for W&M SciClone cluster.
Functions contained here are supposed to be run on SciClone's front-ends
(and perhaps compute nodes as well).
Please see the specific notes in each function.
"""

import os
import os.path
import subprocess
import sys
import tempfile

from wpylib.params.params_struct import Struct as struct
import wpylib.shell_tools as sh

from pyqmc.sites.sciclone import getusername, gethomedir
import pyqmc.sites.sciclone as _site
from pyqmc.utils import _file_search

# Rather general-purpose tools for SciClone:

def fetch_gafqmc_run_output(rundirs, force_update=False, save_walkers=False):
  """Fetches the GAFQMC run output (INFO files and some additional files)
  in bulk.
  Putting them the corresponding result subdirectory for archival.
  Usually this tool is used for runs that fail to finish in time.
  Example usage for Ca+4H2 system:

      >>> fetch_qmc_run_output(glob.glob("part*/rundir"))

  where each rundir is a softlink to the local scratch where the QMC run
  output is stored temporarily.
  """
  from wpylib.regexps import regex
  #destbasedir = tempfile.mkdtemp(dir=os.getcwd())
  #print "Snapshot data is located in subdir: ", destbasedir
  # For now we support only global scratch dir.
  rundir_rx = regex(r'/sciclone/(scr[0-9]+)/' + getusername() + r'/([0-9]+[.-_a-zA-Z0-9]+)\.([-_a-zA-Z0-9]+)\.run')
  if isinstance(rundirs, basestring):
    rundirs = [rundirs]
  files_to_fetch = [
    'INFO',
    'stdout',
    'fort.17',
    'fort.15',
    '*.in',
  ]
  if save_walkers:
    files_to_fetch += [ 'W_old_run_000' ]
  xfiles_to_fetch = ",".join(files_to_fetch)
  for (i, r) in enumerate(rundirs):
    destdir = os.path.dirname(r)
    if destdir == '': destdir = '.'
    destdir_abs = os.path.realpath(destdir)
    rlink = os.path.realpath(r)
    if not (rundir_rx % rlink):
      raise ValueError, "Unparseable dirname: %s" % rlink
    jobid = rundir_rx[2]
    host = rundir_rx[3]
    if not force_update:
      # Check first if the INFO file have the same timestamp:
      old_INFO = _file_search(destdir_abs, ['INFO', 'INFO.gz', 'INFO.bz2', 'INFO.lzma', 'INFO.xz'])
      if old_INFO:
        new_INFO = os.path.join(rlink, 'INFO')
        if os.stat(new_INFO).st_mtime < os.stat(old_INFO).st_mtime:
          print "Skipping: %s -> %s (result already up-to-date)" % (r, rlink)
          continue
    print "Fetching: %s -> %s" % (r, rlink)
    # the rsync flags follow those in run-gafqmc.sh
    sh.run('rsync', ('-ptvb', '%s/{%s}' % (rlink, xfiles_to_fetch), '%s/' % destdir))



# Utilities below relocate stuff from home directory to BigData storage or
# a scratch directory.
# This is for space-saving plus preserving original file paths.

bigfiles_dest_root = {
# Formula for destination root dir for large files (e.g. for
# move_to_* functions below)
  'bigdata': lambda **params: _site.BIGDATA_ROOT + "/%(USER_SCICLONE)s" % params,
  1:         lambda **params: _site.SCR1_ROOT + "/%(USER_SCICLONE)s/BIGFILES" % params,
  2:         lambda **params: _site.SCR2_ROOT + "/%(USER_SCICLONE)s/BIGFILES" % params,
  'trash':   lambda **params: _site.LSCR_ROOT + "/%(USER_SCICLONE)s/deleted" % params,
}

def move_to_bigdata(f):
  """Moves a file (or set of files) from my home dir to BIG_DATA location.
  """
  from os.path import basename, dirname, abspath, realpath, join, exists, islink
  from wpylib.sugar import is_iterable
  if not is_iterable(f):
    f = [f]
  USER_SCICLONE = getusername()
  home_dir = gethomedir()
  lustre_dir = bigfiles_dest_root['lustre'](USER_SCICLONE=USER_SCICLONE)
  for F1 in f:
    if islink(F1): continue
    F1p = realpath(F1)
    F1d = dirname(F1p)
    F1f = basename(F1p)
    if not F1d.startswith(home_dir):
      raise ValueError, \
        "move_to_lustre: path (" + F1p + ") must begin with user's home directory."
    F1d_subpath = F1d[len(home_dir)+1:]
    lustre_path = join(lustre_dir, F1d_subpath)
    #print F1
    #print F1p
    #print lustre_path
    sh.mkdir("-p", "-v", lustre_path)
    sh.provide_link(join(F1d, ".scr"), lustre_path)
    sh.mv("-v", "-i", F1p, join(F1d, ".scr", F1f))
    sh.provide_link(join(F1d, F1f), join(".scr", F1f))


def move_to_scr(f, scr=1):
  """Moves a file (or set of files) from my home dir to /scr1 or /scr2 on PMN.
  """
  assert scr == 1 or scr == 2
  from os.path import basename, dirname, abspath, realpath, join, exists, islink
  from wpylib.sugar import is_iterable
  if not is_iterable(f):
    f = [f]
  USER_SCICLONE = getusername()
  home_dir = gethomedir()
  scr_dir = bigfiles_dest_root[scr](USER_SCICLONE=USER_SCICLONE)
  dotscr = ".scr%d" % scr
  for F1 in f:
    if islink(F1): continue
    F1p = realpath(F1)
    F1d = dirname(F1p)
    F1f = basename(F1p)
    if not F1d.startswith(home_dir):
      raise ValueError, \
        "move_to_scr: path (" + F1p + ") must begin with user's home directory."
    F1d_subpath = F1d[len(home_dir)+1:]
    scr_path = join(scr_dir, F1d_subpath)
    #print F1
    #print F1p
    #print scr_path
    sh.mkdir("-p", "-v", scr_path)
    sh.provide_link(join(F1d, dotscr), scr_path)
    sh.mv("-v", "-i", F1p, join(F1d, dotscr, F1f))
    sh.provide_link(join(F1d, F1f), join(dotscr, F1f))


def move_to_trash(f):
  """Trash a file (or set of files) from my home dir to a designated place
  """
  from os.path import basename, dirname, abspath, realpath, join, exists, islink, isdir
  from wpylib.sugar import is_iterable
  if not is_iterable(f):
    f = [f]
  USER_SCICLONE = getusername()
  home_dir = gethomedir()
  trash_dir = bigfiles_dest_root['trash'](USER_SCICLONE=USER_SCICLONE)
  trash_link = ".trash"
  for F1 in f:
    if islink(F1): continue
    F1p = realpath(F1)
    F1d = dirname(F1p)
    F1f = basename(F1p)
    if not F1d.startswith(home_dir):
      raise ValueError, \
        "move_to_trash: path (" + F1p + ") must begin with user's home directory."
    F1d_subpath = F1d[len(home_dir)+1:]
    trash_path = join(trash_dir, F1d_subpath)
    #print F1
    #print F1p
    #print trash_path
    if not isdir(trash_link):
      if exists(trash_link):
        raise RuntimeError, \
          "Trash link `"+trash_link+"' exists but not a directory in "+F1d
      else:
        sh.mkdir("-p", "-v", trash_path)
        sh.provide_link(join(F1d, trash_link), trash_path)
    sh.mv("-v", "-i", F1p, join(F1d, trash_link, F1f))
    # IMPORTANT: do NOT provide backlink to the trashed file, but only
    # note the subdirectory location:
    #sh.provide_link(join(F1d, F1f), join(".scr", F1f))
