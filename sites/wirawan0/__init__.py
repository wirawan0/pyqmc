#!/usr/bin/python
# $Id: __init__.py,v 1.1 2011-09-12 17:07:29 wirawan Exp $
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

