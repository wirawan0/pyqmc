#!/usr/bin/python
# $Id: gafqmc_quick_dirty.py,v 1.1 2011-09-01 13:36:37 wirawan Exp $
#
# pyqmc.utils.gafqmc_quick_dirty module
#
# Wirawan Purwanto
# Created: 20110831
#

"""
pyqmc.utils.gafqmc_quick_dirty

Rather quick-n-dirty utilities for GAFQMC files and jobs.
"""

_sites_supported = [
  'pmn',
]

import datetime
import numpy
import os
import os.path
import sys
import time

import wpylib.shell_tools as sh
from wpylib.iofmt.text_output import text_output
from wpylib.text_tools import str_grep
from wpylib.sugar import ifelse


def is_qmc_output(filename):
  snippet = sh.pipe_out(("head", "-n", "400", filename), split=True)
  if str_grep("GAFQMC - Generic auxiliary-field quantum Monte Carlo", snippet):
    return True
  elif str_grep("Generic Auxiliary field Quantum Monte Carlo (GAFQMC)", snippet):
    # gen76 and gen79 has this
    return True
  else:
    return False

def is_qmc_output_finished(filename):
  if is_qmc_output(filename):
    snippet = sh.pipe_out(("tail", "-n", "400", filename), split=True)
    if str_grep("Summary of energies:", snippet):
      return True
  return False



def check_qmc(fname=None, archive=False, x=(), xflags=(), reblk=(2,20,2)):
  """Checks the QMC results (INFO file) using my stand-alone QMC inspect tools.
  x is the extra flags to be passed on to check_* tools."""

  if fname == None:
    if os.path.isfile("INFO"):
      fname = ("INFO",)
    else:
      fname = tuple(sh.sorted_glob("part[0-9]*/INFO"))
  elif not hasattr(fname, "__iter__"):
    fname = (fname,)
  else:
    fname = tuple(fname)

  for f in fname:
    if not is_qmc_output(f):
      raise ValueError, "Not a QMC output file: " + f

  if archive:
    farch = open("analysis.txt", "a")
    def prn(x):
      sys.stdout.write(str(x) + "\n")
      farch.write(str(x) + "\n")
  else:
    def prn(x):
      sys.stdout.write(str(x) + "\n")

  prn("")
  prn("Date: " + time.strftime("%Y%m%d"))
  prn("Checking QMC result: " + "\n  + ".join(
          ifelse(len(fname) > 0, [""],[]) + \
          [ os.path.abspath(f) for f in fname ]
      ))
  x = tuple(x) + tuple(xflags)
  if len(x) > 0:
    prn("Extra flags: " + " ".join(x))
  for f in fname:
    if not is_qmc_output_finished(f):
      prn("Warning: calculation unfinished yet: %s" % f)
  # Use the formal tool names here:
  #sh.run("check_eqlb", (fname, "1"))
  prn("EQUILIBRATION/QMC RAW ENERGIES:")
  prn(sh.pipe_out(("check_eqlb",) + fname + ("1",) + x).replace('\x0c', ''))
  prn("REBLOCKING:")
  reblk_info = tuple([ str(vv) for vv in reblk ])
  prn(sh.pipe_out(("check_reblocking",) + fname + reblk_info + ("plt", "txt") + x).replace('\x0c', ''))
  if archive: farch.close()

