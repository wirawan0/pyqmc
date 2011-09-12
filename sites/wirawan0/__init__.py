#!/usr/bin/python
# $Id: __init__.py,v 1.2 2011-09-12 20:47:39 wirawan Exp $
#
# pyqmc.sites.wirawan0 module
#
# Wirawan Purwanto
# Created: 20110912
#

"""
pyqmc.sites.wirawan0

Site-specific support for pyqmc library.
This module is part of PyQMC project.

This module provides support and specific tools for my own computers.

Rather general-purpose tools for this host are available in
pyqmc.sites.wirawan0.tools submodule.

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
site_code = 'wirawan0'

def _detect_site():
  if os.path.isfile("/etc/.wirawan0"):
    return True
  else:
    return False

class site_config(pyqmc.sites.site_config_base):
  """wirawan0-specific configuration."""
  site_code = site_code
  def init_gafqmc(self):
    self._defattr(GAFQMC = os.environ['GAFQMC'])
    self._defattr( \
      GAFQMC_CALC_ROOT = os.path.join(self.GAFQMC, 'qmc'),
    )
    super(site_config, self).init_gafqmc()
  def init_pwqmc(self):
    self._defattr(PWQMC77 = os.environ['PWQMC77'])
    self._defattr( \
      PWQMC_CALC_ROOT = os.path.join(self.PWQMC77, 'qmc'),
    )
    super(site_config, self).init_pwqmc()
