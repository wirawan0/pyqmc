# $Id: programs.py,v 1.1 2011-09-28 19:38:22 wirawan Exp $
#
# pyqmc.nwchem.programs module
#
# Wirawan Purwanto
# Created: 20101028
#

"""
This module contains a general encapsulation of nwchem program.

"""

import sys
import os
import os.path
import time

import wpylib.shell_tools as sh
from wpylib.params import flat as params
from pyqmc.nwchem.output import nwchem_output


class nwchem_program(object):
  """Encapsulation of nwchem program."""
  NWCHEM_EXECUTABLE = "nwchem.sh"
  def __init__(self):
    self.parm = params(
      verbose=2,
      rerun=0,  # forces running even if the output file exists
      parse_results=1,
      logfile=sys.stdout,
    )
  def run(self, infile, **_opts_):
    """Runs nwchem (if the output file hasn't existed) and
    returns the results.

    * verbose (-1): verbosity level
      -  > 2 = even prints out all nwchem output lines to stdout
      -    2 = quiet, but prints the total SCF/DFT energy
      -    1 = prints nothing except saying `nwchem inputfile.nw' to logfile
      - <= 0 = really really quiet

    """
    p = self.parm._create_()
    opt_quiet = []
    nw_exec = self.NWCHEM_EXECUTABLE
    nw_in = infile
    stdout = p.logfile
    verbose = p.verbose

    if nw_in.endswith(".nw"):
      nw_out = nw_in[:-3] + ".out"
    else:
      nw_out = nw_in + ".out"

    if verbose <= 2:
      opt_quiet = ['-q']

    if verbose == 1:
      stdout.write("%s %s -q\n" % (nw_exec,infile,))
      stdout.flush()
    elif verbose == 2 and p.parse_results:
      stdout.write(infile + ":")
      stdout.flush()

    if not os.path.isfile(nw_out):
      sh.run(nw_exec, [nw_in] + opt_quiet)

    if p.parse_results:
      rec = nwchem_output(nw_out)
      if verbose >= 2:
        for e in ('E_SCF', 'E_DFT', 'E_MP2', 'E_CCSD', 'E_CCSD_T',):
          if e in rec:
            stdout.write("%s = %s\n" % (e, rec[e]))
        stdout.flush()
      return rec
