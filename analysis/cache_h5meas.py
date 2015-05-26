#!/usr/bin/python
#
# pyqmc.analysis.cache_h5meas module
#
# Wirawan Purwanto
# Created: 20150522
#

"""
pyqmc.analysis.cache_h5meas
Provides script-wide global cache for measurement HDF5 files.

(WARNING: This is an ad-hoc facility.)
"""

import os
import os.path


# Generalized database for HDF5 measurement objects
# FIXME: This is AD-HOC!!

def H5DB_init():
  global H5DB
  if not "H5DB" in globals():
    H5DB = dict()

def H5DB_close_all():
  global H5DB
  if "H5DB" in globals():
    for v in H5DB.values():
      v.close()
    H5DB.clear()

def H5DB_get(db):
  from pyqmc.results.gafqmc_meas import meas_hdf5 as gmeas_hdf5
  from pyqmc.results.pwqmc_meas import meas_hdf5 as pwmeas_hdf5
  global H5DB
  H5DB_init()
  if isinstance(db, basestring):
    fn_db = os.path.realpath(db)
    if not fn_db in H5DB:
      # FIXME: Use
      H5DB[fn_db] = gmeas_hdf5(fn_db,"r")
    return H5DB[fn_db]
  else:
    return db

