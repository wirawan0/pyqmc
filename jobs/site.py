# -*- python -*-
# $Id: site.py,v 1.1 2011-10-06 19:41:50 wirawan Exp $
#
# pyqmc.jobs.site module
#
# General site-related stuff
#
# Wirawan Purwanto
# Created: 20090320
#

import os

# Root dir, etc. By default they are taken from environment
# variables
# These variables will be None if the variable does not exist
PWQMC77 = os.getenv("PWQMC77")
GAFQMC = os.getenv("GAFQMC")
