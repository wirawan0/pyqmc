# $Id: gafqmc_meas.py,v 1.1 2010-07-29 16:04:19 wirawan Exp $
#
# gafqmc_meas.py
# Tools related to GAFQMC measurement dump (gaf-*.ene) files
#
# Wirawan Purwanto
# Created: 20100608
#
# Based on pwqmc_meas, and using meas_hdf5 and meas_text from that module
# as the base class here.
#

# Python standard modules
import math
import os
import sys
import time

import h5py # HDF5 Python API
import numpy

import pyqmc.results.pwqmc_meas
import pyqmc.results.pwqmc_meas as pwmeas



class meas_hdf5(pwmeas.meas_hdf5):
  """The main object for storing GAFQMC measurement results in HDF5 format.
  """
  header_meta = {
    "Creator": "pyqmc.results.gafqmc_meas",
    "FormatType": "gafqmc_meas",
    "FormatVersion": 0,
    "Devel_CVS_ID": "$Id: gafqmc_meas.py,v 1.1 2010-07-29 16:04:19 wirawan Exp $",
  }


class meas_text(pwmeas.meas_text):
  """A class for reading the measurement dump produced by GAFQMC code
  (i.e. a single pwaf-*.ene file per process).
  The raw dump to be read is in plain text format.
  """
  def check_file_header(self, F):
    L = F.readline().strip()
    if not L.startswith("GAFQMC measurement data: total energies"):
      raise ValueError, "File `" + F.name + "': Invalid file type"


def convert_meas_to_hdf5(output, H0=0, files=None, **opts):
  """A master routine to convert measurement across many gaf-*.meas files to
  the new hdf5 format.
  See documentation of pyqmc.results.pwqmc_meas.convert_meas_to_hdf5
  for more details.

  NOTE: This routine is deprecated.
  Please consider using the more efficient convert_meas_to_hdf5_v2 routine.
  """
  opts.setdefault('meas_hdf5_class', meas_hdf5)
  opts.setdefault('meas_text_class', meas_text)
  return pwmeas.convert_meas_to_hdf5(output, H0=H0, files=files, **opts)


def convert_meas_to_hdf5_v2(output, H0=0, files=None, **opts):
  """A master routine to convert measurement across many gaf-*.meas files to
  the new hdf5 format.
  See documentation of pyqmc.results.pwqmc_meas.convert_meas_to_hdf5_v2
  for more details.
  """
  opts.setdefault('meas_hdf5_class', meas_hdf5)
  opts.setdefault('meas_text_class', meas_text)
  return pwmeas.convert_meas_to_hdf5_v2(output, H0=H0, files=files, **opts)
