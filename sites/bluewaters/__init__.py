#
# pyqmc.sites.bluewaters module
#
# Wirawan Purwanto
# Created: 20130204
# Based on jaguarpf
#

"""
pyqmc.sites.bluewaters

Site-specific support for pyqmc library.
This module is part of PyQMC project.

This module provides support and specific tools for NCSA Blue Waters
supercomputer.

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
site_code = 'bluewaters'
_site_hostnames_sha1 = [
  'e0ece609b72c44e9f7bbd99764c72342c7692472',
  'bc98c167e83039d8e9a9d75d2ad918e2d0bac125',
  '317409ade05086250cee6d74b3d7b4066d6a593c',
  '674cced1ac42a7de1db0ba28fdd8391aad792790'
]
# Note: we use more than one names to check because the names were too
# simple, so we don't want to get fooled.


def _detect_site():
  if pyqmc.sites._hostgrep(_site_hostnames_sha1):
    return True
  else:
    return False

class site_config(pyqmc.sites.site_config_base):
  """Blue Waters-specific configuration."""
  site_code = site_code
  def init(self):
    super(site_config, self).init()
    self.init_bluewaters_dirs()
  def init_bluewaters_dirs(self):
    """Initializes the default BlueWaters-specific dirs:
    - scratch
    - lustre shared space
    """
    pathcat = os.path.join
    self._defattr('SCR_ROOT', pathcat('/scratch/sciteam', self.USER))

def getusername():
  return os.environ['USER']
