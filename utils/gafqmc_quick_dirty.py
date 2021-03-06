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
from wpylib.iofmt.text_input import text_input, head, tail
from wpylib.iofmt.text_output import text_output
from wpylib.text_tools import str_grep
from wpylib.sugar import ifelse


def is_qmc_output(filename):
  snippet = head(filename, 400)
  if str_grep("GAFQMC - Generic auxiliary-field quantum Monte Carlo", snippet):
    return True
  elif str_grep("Generic Auxiliary field Quantum Monte Carlo (GAFQMC)", snippet):
    # gen76 and gen79 has this
    return True
  else:
    return False

def is_qmc_output_finished(filename):
  if is_qmc_output(filename):
    snippet = tail(filename, 400)
    if str_grep("Summary of energies:", snippet):
      return True
  return False



def check_qmc(fname=None, archive=False, xflags=(), eqlb=(1,), reblk=(2,20,2),
              eqlb_flags=(), reblk_flags=(), quiet_stdout=0, force_raw=0):
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

  fname_format = {}
  for f in fname:
    if not is_qmc_output(f):
      if not force_raw:
        raise ValueError, "Not a QMC output file: " + f
      else:
        print >> sys.stderr, "Warning: assuming raw output: " + f
        fname_format[f] = "raw"
    else:
      fname_format[f] = "qmc"

  if archive:
    fnarch = ifelse(isinstance(archive, basestring), archive, "analysis.txt")
    farch = open(fnarch, "a")
    if quiet_stdout:
      def prn(x):
        farch.write(str(x) + "\n")
    else:
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
  for f in fname:
    if fname_format[f] == "qmc" and not is_qmc_output_finished(f):
      prn("Warning: calculation unfinished yet: %s" % f)
  prn("")

  # Use the formal tool names here:
  #sh.run("check_eqlb", (fname, "1"))
  if eqlb:
    prn("EQUILIBRATION/QMC RAW ENERGIES:")
    eqlb_info = tuple(map(str, eqlb))
    x = tuple(xflags) + tuple(eqlb_flags)
    if True or len(x) > 0:
      prn("Flags: " + " ".join(eqlb_info + x))
    prn(sh.pipe_out(("check_eqlb",) + fname + eqlb_info + x).replace('\x0c', ''))
  
  if reblk:
    prn("REBLOCKING:")
    reblk_info = tuple(map(str, reblk))
    x = tuple(xflags) + tuple(reblk_flags)
    if True or len(x) > 0:
      prn("Flags: " + " ".join(reblk_info + x))
    prn(sh.pipe_out(("check_reblocking",) + fname + reblk_info + ("plt", "txt") + x).replace('\x0c', ''))
  if archive: farch.close()

