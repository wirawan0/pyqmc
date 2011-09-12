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
_site_hostname_sha1 = 'f1e7e1de72c137f5e286567137f2c9af5fc92964'


def _detect_site():
  if str_grep(_site_hostname_sha1, pyqmc.sites._etc_hosts_ipv4hosts_digest):
    return True
  else:
    return False

class site_config(pyqmc.sites.site_config_base):
  """PMN-specific configuration."""
  site_code = site_code
