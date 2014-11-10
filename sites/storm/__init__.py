#!/usr/bin/python
#
# pyqmc.sites.storm module
#
# Wirawan Purwanto
# Created: 20141030
#

"""
pyqmc.sites.storm

Site-specific support for pyqmc library.
This module is part of PyQMC project.

This module provides support and specific tools for W&M's storm cluster.

Rather general-purpose tools for storm are available in
pyqmc.sites.storm.tools submodule.
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
site_code = 'storm'
_site_hostname_sha1   = [
  '8927a74a85d53840f6e7986d3c200b2239915413',
]
_compute_node_prefixes = [
  'ra', 'sn', 'ha', 'ice', 'wi',
]

LSCR_ROOT = '/local/scr'
SCR1_ROOT = '/sciclone/scr00'
SCR2_ROOT = '/sciclone/scr1'  # TENTATIVE
# other scr partitions exist, see sciclone manual page
BIGDATA_ROOT = '/storm/data10' # the former /mnt/lustre

def _detect_site():
  from pyqmc.sites import _etc_hosts_ipv4hosts_digest, _info
  for s in _site_hostname_sha1:
    if str_grep(s, _etc_hosts_ipv4hosts_digest):
      # This alone used to be a good indicator of sciclone,
      # but not any more since the introduction of storm etc
      # as part of W&M HPC system.
      #
      # So we need to add several further checks:
      # - Front-end check: against the FQDN
      if _info.fqdn_sha1 in _site_hostname_sha1:
        return True
      # - Compute-node check: against the literal short name
      #   (FIXME: This is rather weak check; maybe useful to think of another way)
      elif _info.hostshortname.rstrip("0123456789") in _compute_node_prefixes:
        return True
  return False

class site_config(pyqmc.sites.site_config_base):
  """Storm-specific configuration."""
  site_code = site_code
  def init(self):
    super(site_config, self).init()
    self.init_storm_dirs()
  def init_storm_dirs(self):
    """Initializes the default Storm-specific dirs:
    - scratch space
    - shared space
    """
    pathcat = os.path.join
    self._defattr('LSCR_ROOT', pathcat(LSCR_ROOT, self.USER))
    self._defattr('SCR1_ROOT', pathcat(SCR1_ROOT, self.USER))
    self._defattr('SCR2_ROOT', pathcat(SCR2_ROOT, self.USER))
    self._defattr('BIGDATA_ROOT', pathcat(BIGDATA_ROOT, self.USER))
  def init_gafqmc(self):
    pathcat = os.path.join
    self.init_storm_dirs()
    self._defattr( \
      GAFQMC_BIG_ROOT = pathcat(self.BIGDATA_ROOT, "BIGFILES.runtime/GAFQMC"),
      NWCHEM_BIG_ROOT = pathcat(self.BIGDATA_ROOT, "BIGFILES.runtime/GAFQMC/nwchem"),
      GAMESS_BIG_ROOT = pathcat(self.BIGDATA_ROOT, "BIGFILES.runtime/GAFQMC/gamess"),
      GAUSSIAN_BIG_ROOT = pathcat(self.BIGDATA_ROOT, "BIGFILES.runtime/GAFQMC/gaussian"),
    )
    super(site_config, self).init_gafqmc()
  def init_pwqmc(self):
    pathcat = os.path.join
    self.init_storm_dirs()
    self._defattr( \
      PWQMC_BIG_ROOT = pathcat(self.BIGDATA_ROOT, "BIGFILES.runtime/PWQMC-77"),
    )
    super(site_config, self).init_pwqmc()

def getusername():
  """Override this function if needed (e.g. not running on the site)."""
  return os.environ['USER']

def gethomedir():
  """Override this function if needed (e.g. not running on the site)
  to return the user's home directory in sciclone."""
  return os.environ['HOME']
