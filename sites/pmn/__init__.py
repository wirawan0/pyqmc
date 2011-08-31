#!/usr/bin/python
# $Id: __init__.py,v 1.1 2011-08-31 20:15:26 wirawan Exp $
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

import wpylib.shell_tools as sh


# Standard information variables
site_code = 'pmn'


