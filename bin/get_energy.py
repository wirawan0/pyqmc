#!/usr/bin/env python
#
# get_energy.py
#
# Wirawan Purwanto
# Created: 20130220
#

"""
pyqmc.bin.get_energy

Utility to extract energies from perform some operation to it (e.g.
averaging, etc).
"""

import datetime
import numpy
import os
import os.path
import sys
import time

from wpylib.text_tools import str_grep
import wpylib.shell_tools as sh
from wpylib.iofmt.text_input import text_input, head, tail
from pyqmc.results.gafqmc_info import gafqmc_info
from pyqmc.results.pwqmc_info import pwqmc_info
from pyqmc.results.meas_merge import merged_measurement


def is_afqmc_output(filename):
  snippet = head(filename, 400)
  if str_grep("GAFQMC - Generic auxiliary-field quantum Monte Carlo", snippet):
    return "GAFQMC"
  elif str_grep("Generic Auxiliary field Quantum Monte Carlo (GAFQMC)", snippet):
    # gen76 and gen79 has this
    return "GAFQMC"
  elif str_grep("Planewave-AFQMC calculation for system: General electronic system", snippet):
    # gen76 and gen79 has this
    return "PWQMC"
  else:
    return False


def get_energy_merged(files):
  """An alternative of get-energy.gawk --merged tool.
  """
  recs = []
  for f in files:
    ftype = is_afqmc_output(f)
    if ftype == "GAFQMC":
      info = gafqmc_info(f)
      print "# Reading GAFQMC info file: ", f
    elif ftype == "PWQMC":
      info = pwqmc_info(f)
      print "# Reading PWQMC info file:  ", f
    else:
      raise PyqmcDataError, \
        "Invalid input file type: %s: not a supported AFQMC INFO file" % (f,)
    recs.append(info)
  Merge = merged_measurement()
  Merge.merge([ I.meas_energy for I in recs ])
  Merge.print_formatted("%(beta)5.2f %(overlap)18.9f %(Etotal)15.9f %(ndata)d\n")


if __name__ == "__main__":
  get_energy_merged(sys.argv[1:])
