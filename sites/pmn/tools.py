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
"""

import os
import os.path
import subprocess
import sys
import tempfile

from wpylib.params.params_struct import Struct as struct
import wpylib.shell_tools as sh

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
  rundir_rx = regex(r'/state/partition1/' + os.environ['USER'] + '/([0-9]+)\.([-_a-zA-Z0-9]+)\.run')
  if isinstance(rundirs, str):
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
  rundir_rx = regex(r'/state/partition1/' + os.environ['USER'] + '/([0-9]+)\.([-_a-zA-Z0-9]+)\.run')
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

