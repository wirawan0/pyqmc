#!/usr/bin/python
# $Id: __init__.py,v 1.3 2011-09-12 20:47:37 wirawan Exp $
#
# pyqmc.sites.avocado module
#
# Wirawan Purwanto
# Created: 20110831
#

"""
pyqmc.sites.avocado

Site-specific support for pyqmc library.
This module is part of PyQMC project.

This module provides support and specific tools for a W&M condensed
matter group's private cluster.

Rather general-purpose tools for AVOCADO are available in
pyqmc.sites.avocado.tools submodule.

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
site_code = 'avocado'
_site_hostname_sha1 = '6f9ae820bce771e46a6da011820afce889c87416'


def _detect_site():
  if str_grep(_site_hostname_sha1, pyqmc.sites._etc_hosts_ipv4hosts_digest):
    return True
  else:
    return False

class site_config(pyqmc.sites.site_config_base):
  """AVOCADO-specific configuration."""
  site_code = site_code
  def init(self):
    super(site_config, self).init()
    self.init_avocado_dirs()
  def init_avocado_dirs(self):
    """Initializes the default AVOCADO-specific dirs:
    - scratch
    - lustre shared space
    """
    pathcat = os.path.join
    self._defattr('LSCR_ROOT', pathcat('/lscr', self.USER))

