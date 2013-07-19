#
# pyqmc.sites.jaguarpf module
#
# Wirawan Purwanto
# Created: 20120613
#

"""
pyqmc.sites.jaguarpf

Site-specific support for pyqmc library.
This module is part of PyQMC project.

This module provides support and specific tools for OLCF jaguarpf
supercomputer and its predecessors (titan, ...).

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
site_code = 'jaguarpf'
_site_hostnames_sha1 = [
  'aa9c409cfbc3d66e77bdcf5d8f6693b923bec76f',
  '532280b8b4297b4869c7d82f46f9c310f536a433',
]


def _detect_site():
  if pyqmc.sites._hostgrep(_site_hostnames_sha1):
    return True
  else:
    return False

class site_config(pyqmc.sites.site_config_base):
  """JaguarPF-specific configuration."""
  site_code = site_code
  def init(self):
    super(site_config, self).init()
    self.init_jaguarpf_dirs()
  def init_jaguarpf_dirs(self):
    """Initializes the default JaguarPF-specific dirs:
    - scratch
    - lustre shared space
    """
    pathcat = os.path.join
    self._defattr('SCR_ROOT', pathcat('/tmp/work', self.USER))

def getusername():
  return os.environ['USER']
