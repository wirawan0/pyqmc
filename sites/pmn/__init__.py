#!/usr/bin/python
# $Id: __init__.py,v 1.2 2011-09-01 20:10:37 wirawan Exp $
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
  if strgrep(_site_hostname_sha1, pyqmc.sites._etc_hosts_ipv4hosts_digest):
    return True
  else:
    return False

