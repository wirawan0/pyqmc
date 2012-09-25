#!/usr/bin/python
# $Id: tools.py,v 1.1 2011-08-31 20:15:26 wirawan Exp $
#
# pyqmc.sites.pmn.tools module
#
# Wirawan Purwanto
# Created: 20110831
#

"""
pyqmc.sites.pmn.tools

Site-specific support for pyqmc library.
This module is part of PyQMC project.

This module provides specific tools for CPD's PMN Opteron cluster.
Functions contained here are supposed to be run on PMN front-end
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

from pyqmc.sites.pmn import getusername

# Rather general-purpose tools for PMN:

def get_gafqmc_run_snapshot(rundirs):
  """Fetches the GAFQMC run snapshot (the INFO files) in bulk.
  Puts them all in temporary subdirectory for subsequent examination.
  This tool is used for runs started by the run-gafqmc.sh script.
  Example usage for Ca+4H2 system:

      >>> get_qmc_run_snapshot(glob.glob("part*/rundir"))

  where each rundir is a softlink to local scratch where the QMC run
  output is stored temporarily.
  """
  from wpylib.regexps import regex
  from os.path import abspath, dirname, join
  destbasedir = tempfile.mkdtemp(dir=os.getcwd())
  print "Snapshot data is located in subdir: ", destbasedir
  # This is the standard location for rundir on local scratch as
  # defined by run-gafqmc.sh:
  rundir_rx = regex(r'/state/partition1/' + getusername() + '/([0-9]+)\.([-_a-zA-Z0-9]+)\.run')
  if isinstance(rundirs, basestring):
    rundirs = [rundirs]
  for (i, r) in enumerate(rundirs):
    destdir = os.path.join(destbasedir, "part%04d" % i)
    os.mkdir(destdir)
    os.symlink(abspath(dirname(r)), join(destdir,"rsltdir"))
    rlink = os.path.realpath(r)
    if not (rundir_rx % rlink):
      raise ValueError, "Unparseable dirname: %s" % rlink
    jobid = rundir_rx[1]
    host = rundir_rx[2]
    print "Fetching: %s -> %s" % (r, rlink)
    sh.run('scp', ('-p', '%s:%s/INFO' % (host, rlink), '%s/INFO' % destdir))


def fetch_gafqmc_run_output(rundirs, force_update=False, save_walkers=False):
  """Fetches the GAFQMC run snapshot (INFO files and some additional files)
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
  rundir_rx = regex(r'/state/partition1/' + getusername() + '/([0-9]+)\.([-_a-zA-Z0-9]+)\.run')
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
    jobid = rundir_rx[1]
    host = rundir_rx[2]
    if not force_update:
      # Check first if the INFO file have the same timestamp:
      old_INFO = _file_search(destdir_abs, ['INFO', 'INFO.lzma'])
      if old_INFO:
        new_INFO = os.path.join(rlink, 'INFO')
        if cnode_shell(host, 'test %s -nt %s' % (new_INFO, old_INFO)) != 0:
          print "Skipping: %s -> %s (result already up-to-date)" % (r, rlink)
          continue
    print "Fetching: %s -> %s" % (r, rlink)
    # the rsync flags follow those in run-gafqmc.sh
    sh.run('rsync', ('-ptvb', '%s:%s/{%s}' % (host, rlink, xfiles_to_fetch), '%s/' % destdir))


def cnode_shell(host, args, **opts):
  """Invokes a command on compute node, returning the error code.
  Beware: the arguments must be a string, and it *will* be shell-expanded
  upon execution on the remote node.
  """
  return subprocess.call(('ssh', host, args), **opts)


def _file_search(Dir, Names):
  for N in Names:
    fn = os.path.join(Dir, N)
    if os.path.isfile(fn):
      return fn
  return None


def get_running_jobs(sort=True):
  """Obtains the list of currently running jobs for the user."""
  username = getusername()
  # see squeue man page for status code (%t specifier)
  listjob = sh.pipe_out(("squeue", "-u", username, "--noheader", "--format=%i %t"), split=True)
  rslt = []
  # was using `badstatus' this on my awk scriplet, but I changed to `runstatus'
  # below for positive status checking.
  badstatus = (
    'PD', # pending
    'CA', # cancelled
    'CD', # completed
  )
  # treat one of these statuses as "running"
  runstatus = (
    'CF', # configuring
    'R',  # running
    'F',  # failure
    'NF', # node failure
    'TO', # timeout
    'CG', # completing
  )
  for job1 in listjob:
    R = job1.split()
    if R[1] in runstatus:
      rslt.append(R[0])
  if sort:
    rslt.sort()
  return rslt


def get_pending_jobs(sort=True):
  """Obtains the list of currently pending (queued) jobs for the user."""
  username = getusername()
  # see squeue man page for status code (%t specifier)
  listjob = sh.pipe_out(("squeue", "-u", username, "--noheader", "--format=%i %t"), split=True)
  rslt = []
  # treat one of these statuses as "running"
  qstatus = (
    'PD', # pending
    'S',  # suspended
  )
  for job1 in listjob:
    R = job1.split()
    if R[1] in qstatus:
      rslt.append(R[0])
  if sort:
    rslt.sort()
  return rslt


# Utilities below relocate stuff from home directory to lustre or
# a scratch directory.
# This is for space-saving plus preserving original file paths.

bigfiles_dest_root = {
# Formula for destination root dir for large files (e.g. for
# move_to_* functions below)
  'lustre': lambda **params: "/mnt/lustre/%(USER_PMN)s/BIGFILES.runtime" % params,
  1:        lambda **params: "/scr1/%(USER_PMN)s/BIGFILES" % params,
  2:        lambda **params: "/scr2/%(USER_PMN)s/BIGFILES" % params,
  'trash':  lambda **params: "/scr1/%(USER_PMN)s/deleted" % params,
}

def move_to_lustre(f):
  """Moves a file (or set of files) from my home dir to lustre PFS on PMN.
  """
  from os.path import basename, dirname, abspath, realpath, join, exists, islink
  from wpylib.sugar import is_iterable
  if not is_iterable(f):
    f = [f]
  USER_PMN = getusername()
  home_dir = join("/home", USER_PMN)
  lustre_dir = bigfiles_dest_root['lustre'](USER_PMN=USER_PMN)
  # join("/mnt/lustre", USER_PMN, "BIGFILES.runtime")
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


def move_to_scr(f, scr=2):
  """Moves a file (or set of files) from my home dir to /scr1 or /scr2 on PMN.
  """
  assert scr == 1 or scr == 2
  from os.path import basename, dirname, abspath, realpath, join, exists, islink
  from wpylib.sugar import is_iterable
  if not is_iterable(f):
    f = [f]
  USER_PMN = getusername()
  home_dir = join("/home", USER_PMN)
  scr_dir = bigfiles_dest_root[scr](USER_PMN=USER_PMN)
  # join("/scr%d" % scr, USER_PMN, "BIGFILES")
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
  USER_PMN = getusername()
  home_dir = join("/home", USER_PMN)
  trash_dir = bigfiles_dest_root['trash'](USER_PMN=USER_PMN)
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
