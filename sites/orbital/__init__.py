#!/usr/bin/python
#
# pyqmc.sites.orbital module
#
# Wirawan Purwanto
# Created: 20140224
#

"""
pyqmc.sites.orbital

Site-specific support for pyqmc library.
This module is part of PyQMC project.

This module provides support and specific tools for my own computers.

Rather general-purpose tools for this host are available in
pyqmc.sites.orbital.tools submodule.

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
site_code = 'orbital'
_site_hostname_sha1 = '1e0e7e9b3d26ba8663eb9c64692383dd1c29eef9'

def _detect_site():
  if str_grep(_site_hostname_sha1, pyqmc.sites._etc_hosts_ipv4hosts_digest):
    return True
  else:
    return False

class site_config(pyqmc.sites.site_config_base):
  """orbital-specific configuration."""
  site_code = site_code
  def init_gafqmc(self):
    self._defattr(GAFQMC = os.environ['GAFQMC'])
    self._defattr( \
      GAFQMC_CALC_ROOT = os.path.join(self.GAFQMC),
    )
    super(site_config, self).init_gafqmc()
  def init_pwqmc(self):
    self._defattr(PWQMC77 = os.environ['PWQMC77'])
    self._defattr( \
      PWQMC_CALC_ROOT = os.path.join(self.PWQMC77),
    )
    super(site_config, self).init_pwqmc()
