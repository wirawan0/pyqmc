# $Id: gafqmc_info.py,v 1.3 2011-03-09 15:44:47 wirawan Exp $
#
# gafqmc_info.py
# Tools to parse GAFQMC INFO file
#
# Wirawan Purwanto
# Created: 20101025
#
# IMPORTANT: Try to make field names consistent with those in pwqmc_info.

# Python standard modules
import math
import os
import os.path
import re
import sys
import time

import numpy

from wpylib.iofmt.text_input import text_input
from wpylib.db.result_base import result_base
from wpylib.sugar import ifelse

class gafqmc_info(result_base):
  '''Structure to represent the metadata contained in INFO file
  (GAFQMC version).

  Available information:
  * info_file
  * start_time
  * info_mtime
  * calc_time (defined as info_mtime - start_time in seconds)
  * nbasis
  * Evar_noconst  Evar  H0
  * Etrial_noconst  Etrial
  * deltau  betablk  nblkstep
  * nwlk nwlkmax nwlkmin
  * itv_Em itv_pc itv_pc_eq

  '''
  meas_dtype = numpy.dtype([('beta',float), ('overlap',float),
                            ('Etotal',float), ('Eproj',float)])
  def parse_INFO(self, INFO):
    '''Gets all the necessary info (calculation parameters) from a
    GAFQMC INFO file.
    This is a very old routine.
    We use this as temporary starting point.'''
    from pyqmc import PyqmcParseError

    info_file = text_input(INFO)
    self.clear()
    rslt = self
    rslt['info_file'] = INFO
    rslt['info_mtime'] = time.localtime(os.path.getmtime(INFO))
    for L in info_file:
      Ls = L.strip()
      flds = Ls.split()
      if len(flds) == 0:
        continue
      elif flds[0] == "Variational":
        if flds[1] == "energy":
          rslt["Evar"] = float(flds[3])
        elif flds[1] == "energy=":
          rslt["Evar"] = float(flds[2])
      elif flds[0] == "nbasis":
        rslt["nbasis"] = int(flds[2])
      elif flds[0] == "Energy_N_QMC":
        rslt["H0"] = float(flds[1])
      elif flds[0] == "deltau=":
        rslt["deltau"] = float(flds[1])
      elif flds[0] == "beta=":
        rslt["betablk"] = float(flds[1])
      elif Ls.startswith("input etrial="):
        rslt["Etrial_noconst"] = float(flds[2]) # no H0 yet
        # anorm is also available on the same line:
        rslt["anorm"] = float(flds[5])
      elif Ls.startswith("New etrial to be used in El_bound:"):
        rslt["Etrial_noconst"] = float(flds[8]) # H0 specified below
      elif Ls.startswith("nblk="):
        rslt["nblk"] = int(flds[1])
      elif Ls.startswith("neq="):
        rslt["neq"] = int(flds[1])
      elif Ls.startswith("ngrth="):
        rslt["ngrth"] = int(flds[1])
      elif Ls.startswith("No Growth phase:"):
        rslt["ngrth"] = 0
      elif Ls.startswith("itv_em="):
        rslt["itv_em"] = int(flds[1])
      elif Ls.startswith("itv_pc="):
        rslt["itv_pc"] = int(flds[1])
      elif Ls.startswith("itv_pc_eq="):
        rslt["itv_pc_eq"] = int(flds[1])
      elif Ls.startswith("nblkstep="):
        rslt["nblkstep"] = int(flds[1])
      elif Ls.startswith("nwlk="):
        rslt["nwlk"] = int(flds[1])
      elif Ls.startswith("nwlkmax="):
        rslt["nwlkmax"] = int(flds[1])
      elif Ls.startswith("nwlkmin="):
        rslt["nwlkmin"] = int(flds[1])

      # ---runtime info below---
      elif Ls.startswith("Using OpenMP with number of threads = "):
        rslt["num_threads"] = int(flds[7])
      elif Ls.startswith("Parallel version of GAFQMC, using NProc = "):
        rslt["code_name"] = "gafqmc"
        rslt["code_branch"] = "mpi"
        rslt["num_tasks"] = int(flds[7])
      elif Ls.startswith("Host:"):
        rslt["run_host"] = flds[1]
      elif Ls.startswith("Program was run on"):
        # CAVEAT: start_time can be off by several hour if the local time
        # zone is different from the time zone where the calculation
        # was done.
        # FIXME this!
        rslt['start_time'] = \
          time.strptime(flds[4] + " " + flds[6], "%Y/%m/%d %H:%M:%S")
      elif Ls.startswith("Program was ended on"):
        rslt['end_time'] = \
          time.strptime(flds[4] + " " + flds[6], "%Y/%m/%d %H:%M:%S")

      # measurement and other complex data capture
      elif Ls.startswith("Measurement") and flds[1].startswith("phase...."):
        #print "found meas!"
        self.locate_text_marker(info_file,
          (lambda S : S.startswith("Output:")),
          max_try=30,
          errmsg="Cannot locate the beginning of measurement data")
        self.parse_measurement0(info_file, rslt)

    rslt.setdefault("nwlkmax", rslt.nwlk * 2)
    rslt.setdefault("nwlkmin", max(rslt.nwlk / 2, 1))

    rslt["Evar_noconst"] = rslt["Evar"] - rslt["H0"]
    rslt["Etrial"] = rslt["Etrial_noconst"] + rslt["H0"]
    rslt["calc_time"] = time.mktime(rslt[ifelse("end_time" in rslt, "end_time", "info_mtime")]) \
                      - time.mktime(rslt["start_time"])
    rslt["run_mpi"] = ("num_tasks" in rslt)
    rslt["run_openmp"] = ("num_threads" in rslt)
    return rslt

  def locate_text_marker(self, info_file, match_func, max_try, errmsg):
    """Seeks the text lines until a given marker is found.
    An exception is raised if after max_try read attempts,
    the marker is not found."""
    for i in xrange(max_try):
      Ls = info_file.next().strip()
      if match_func(Ls):
        return True
    raise PyqmcParseError, errmsg

  def parse_measurement0(self, info_file, rslt):
    """Internal routine to parse only the measurement results of the
    file.
    info_file is an open file-like object.
    The last line read must have been 'Measurement phase...'

    TODO:
    - add stand-alone parse_measurement routine?
    """
    from pyqmc import PyqmcParseError

    # FIXME: Add beginning marker detection (previous text line must be
    # "Output:")

    for_D2E = lambda s : s.replace("D","E").replace("d","e")
    EOS = re.compile(r"^\s*Final Results:\s*$")  # end-of-stream marker
    RS = re.compile(r"^\s*-+\s*$")  # record separator

    meas = []
    for L in info_file:
      Ls = for_D2E(L.strip())
      flds = Ls.split()
      if EOS.search(Ls):
        break   # end-of-stream detected
      elif len(flds) == 3:
        # special case to handle wrapped Fortran output
        Ls2 = for_D2E(info_file.next().strip())
        flds2 = Ls2.split()
        if len(flds2) == 0:
          raise PyqmcParseError, \
            "Invalid format in GAFQMC measurement text (INFO)"
        flds.append(flds2[0])
      elif len(flds) < 4:
        raise PyqmcParseError, \
          "Invalid format in GAFQMC measurement text (INFO)"
      try:
        rec = tuple(map((lambda x: float(x.rstrip(','))), flds[:4]))
      except:
        raise PyqmcParseError, \
          "Error parsing GAFQMC measurement text (INFO)"
      meas.append(rec)

      self.locate_text_marker(info_file,
        (lambda S : RS.search(S)), max_try=20,
        errmsg="Cannot locate a valid record separator in GAFQMC measurement text (INFO)")

    dtype = self.meas_dtype
    rslt["meas_energy"] = numpy.array(meas, dtype=dtype)

  parse_file_ = parse_INFO
