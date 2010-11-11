# $Id: gafqmc_info.py,v 1.2 2010-11-11 18:03:46 wirawan Exp $
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
import sys
import time

import numpy

from wpylib.iofmt.text_input import text_input
from wpylib.db.result_base import result_base

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
  def parse_INFO(self, INFO):
    '''Gets all the necessary info (calculation parameters) from a
    GAFQMC INFO file.
    This is a very old routine.
    We use this as temporary starting point.'''
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
      elif Ls.startswith("Program was run on"):
        # CAVEAT: start_time can be off by 1 hour if the local time
        # zone is different from the time zone where the calculation
        # was done.
        # FIXME this!
        rslt['start_time'] = \
          time.strptime(flds[4] + " " + flds[6], "%Y/%m/%d %H:%M:%S")
      elif flds[0] == "Variational" and flds[1] == "energy":
        rslt["Evar"] = float(flds[3])
      elif flds[0] == "nbasis":
        rslt["nbasis"] = int(flds[2])
      elif flds[0] == "Energy_N_QMC":
        rslt["H0"] = float(flds[1])
      elif flds[0] == "deltau=":
        rslt["deltau"] = float(flds[1])
      elif flds[0] == "beta=":
        rslt["betablk"] = float(flds[1])
      elif Ls.startswith("Input etrial="):
        rslt["Etrial_noconst"] = float(flds[2]) # no H0 yet
        # anorm is also available on the same line:
        rslt["anorm"] = float(flds[5])
      elif Ls.startswith("New etrial to be used in El_bound:"):
        rslt["Etrial_noconst"] = float(flds[8]) # H0 specified below
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
    rslt.setdefault("nwlkmax", rslt.nwlk * 2)
    rslt.setdefault("nwlkmin", max(rslt.nwlk / 2, 1))

    rslt["Evar_noconst"] = rslt["Evar"] - rslt["H0"]
    rslt["Etrial"] = rslt["Etrial_noconst"] + rslt["H0"]
    rslt["calc_time"] = time.mktime(rslt["info_mtime"]) \
                      - time.mktime(rslt["start_time"])
    return rslt

  parse_file_ = parse_INFO
