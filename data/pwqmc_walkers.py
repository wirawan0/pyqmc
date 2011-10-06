#!/usr/bin/python
# $Id: pwqmc_walkers.py,v 1.1 2011-10-06 19:41:50 wirawan Exp $
#
# pyqmc.data.pwqmc_walkers module
#
# Tools for handling PWQMC-77 walker files
#
# Wirawan Purwanto
# Created: 20090320
#

class wlk_pwaf(object):
    '''PWAF-style walker files.
    Each MPI process (or head of the process group in the split-walker
    case) has its own pwaf-NNNNN.wlk file.
    '''
    # By default assume the walkers reside in a subdir called "wlk/"
    wlkdir = "wlk"


