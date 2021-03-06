#!/usr/bin/python
# $Id: __init__.py,v 1.3 2011-09-12 20:47:37 wirawan Exp $
#
# pyqmc.sites.pmn module
#
# Wirawan Purwanto
# Created: 20110831
#

"""
pyqmc.sites.pmn

Site-specific support for pyqmc library.
This module is part of PyQMC project.

This module provides support and specific tools for CPD's PMN Opteron
cluster.

Rather general-purpose tools for PMN are available in
pyqmc.sites.pmn.tools submodule.

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
site_code = 'pmn'
_site_hostnames_sha1 = [
  '10fe1b1d2a89fa3587b142b3722c5eef73231a2e',
  'f1e7e1de72c137f5e286567137f2c9af5fc92964',
]


def _detect_site():
  import re
  from pyqmc.sites import _info
  if pyqmc.sites._hostgrep(_site_hostnames_sha1):
      # This alone used to be a good indicator of pmn,
      # but not any more since the introduction of storm
      # as part of W&M HPC system.
      #
      # So we need to add several further checks:
      # - Front-end check: against the FQDN
      if _info.fqdn_sha1 in _site_hostnames_sha1:
        return True
      # - Compute-node check: against the literal domain name: "cNNN.local"
      #   (FIXME: This is rather weak check; it maybe useful to think of another way)
      elif re.search(r'^c[0-9]+\.local$', _info.hostname):
        return True
  else:
    return False

class site_config(pyqmc.sites.site_config_base):
  """PMN-specific configuration."""
  site_code = site_code
  def init(self):
    super(site_config, self).init()
    self.init_pmn_dirs()
  def init_pmn_dirs(self):
    """Initializes the default PMN-specific dirs:
    - scratch
    - lustre shared space
    """
    pathcat = os.path.join
    self._defattr('LSCR_ROOT', pathcat('/lscr', self.USER))
    self._defattr('SCR1_ROOT', pathcat('/scr1', self.USER))
    self._defattr('SCR2_ROOT', pathcat('/scr2', self.USER))
    self._defattr('LUSTRE_ROOT', pathcat('/mnt/lustre', self.USER))
  def init_gafqmc(self):
    pathcat = os.path.join
    self.init_pmn_dirs()
    self._defattr( \
      GAFQMC_BIG_ROOT = pathcat(self.LUSTRE_ROOT, "BIGFILES.runtime/GAFQMC"),
      NWCHEM_BIG_ROOT = pathcat(self.LUSTRE_ROOT, "BIGFILES.runtime/GAFQMC/nwchem"),
      GAMESS_BIG_ROOT = pathcat(self.LUSTRE_ROOT, "BIGFILES.runtime/GAFQMC/gamess"),
      GAUSSIAN_BIG_ROOT = pathcat(self.LUSTRE_ROOT, "BIGFILES.runtime/GAFQMC/gaussian"),
    )
    super(site_config, self).init_gafqmc()
  def init_pwqmc(self):
    pathcat = os.path.join
    self.init_pmn_dirs()
    self._defattr( \
      PWQMC_BIG_ROOT = pathcat(self.LUSTRE_ROOT, "BIGFILES.runtime/PWQMC-77"),
    )
    super(site_config, self).init_pwqmc()

def getusername():
  return os.environ['USER']
