# -*- python -*-
# $Id: acct.py,v 1.1 2011-10-06 19:41:50 wirawan Exp $
#
# pyqmc.jobs.acct module
#
# Job accounting
#
# Wirawan Purwanto
# Created: 20090320
#

import pyqmc.jobs.site
from pyqmc.jobs.site import PWQMC77, GAFQMC


class Acct(object):
    '''Job accounting master object.'''
    def submit(self, jobid, site_info, calc_info):
        '''Submit an entry of the job into the database.'''

