#!/usr/bin/python
#
# pyqmc.sites.sciclone module
#
# Wirawan Purwanto
# Created: 20120925
#

"""
pyqmc.sites.sciclone

Site-specific support for pyqmc library.
This module is part of PyQMC project.

This module provides support and specific tools for W&M's SciClone
cluster.

Rather general-purpose tools for SciClone are available in
pyqmc.sites.sciclone.tools submodule.

"""

import os
import os.path
import sys
import tempfile
import pyqmc.sites

import wpylib.shell_tools as sh

import pyqmc.sites
from pyqmc.sites import str_grep


# Standard information variables
site_code = 'sciclone'
_site_hostname_sha1   = [
  'c2157b99683f8d4feca2440f9df024569334e2ad',
  '6c79982b1c616b68625a61006de37f0e2f1dc6a7',
  '2f22b7c85e1f144eb0005e887236cc4451c20cb7',
]

LSCR_ROOT = '/local/scr'
SCR1_ROOT = '/sciclone/scr20'
SCR2_ROOT = '/sciclone/scr02'
# other scr partitions exist, see sciclone manual page
BIGDATA_ROOT = '/sciclone/data10'

def _detect_site():
  for s in _site_hostname_sha1:
    if str_grep(s, pyqmc.sites._etc_hosts_ipv4hosts_digest):
      return True
  return False

class site_config(pyqmc.sites.site_config_base):
  """SciClone-specific configuration."""
  site_code = site_code
  def init(self):
    super(site_config, self).init()
    self.init_sciclone_dirs()
  def init_sciclone_dirs(self):
    """Initializes the default PMN-specific dirs:
    - scratch space
    - shared space
    """
    pathcat = os.path.join
    self._defattr('LSCR_ROOT', pathcat(LSCR_ROOT, self.USER))
    self._defattr('SCR1_ROOT', pathcat(SCR1_ROOT, self.USER))
    self._defattr('SCR2_ROOT', pathcat(SCR2_ROOT, self.USER))
    self._defattr('BIGDATA_ROOT', pathcat(BIGDATA_ROOT, self.USER))

def getusername():
  """Override this function if needed (e.g. not running on the site)."""
  return os.environ['USER']

def gethomedir():
  """Override this function if needed (e.g. not running on the site)
  to return the user's home directory in sciclone."""
  return os.environ['HOME']
