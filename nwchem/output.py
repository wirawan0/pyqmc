# $Id: output.py,v 1.2 2010-12-01 17:09:39 wirawan Exp $
#
# pyqmc.nwchem.output module
#
# Wirawan Purwanto
# Created: 20101028
#

"""
This module contains an output parser for nwchem code.

"""

import os
import os.path
import time

from wpylib.regexps import regex
from wpylib.iofmt.text_input import text_input
from wpylib.db.result_base import result_base

class nwchem_output(result_base):
  """Parser and structure for nwchem output text file.

  Fields from results_base:
    * filename_
    * absfilename_ -- full file name including absolute directory path.
    Here:
  Fields defined here:
    * info_code_version  nwchem version used in the calculation
    * info_mtime
  """

  def parse_file_(self, filename):
    """Extracts information from an nwchem output file (from its stdout).
    Right now, this parser is only good for single-point calculations
    (i.e. no multijob or geometry optimization at this point)."""
  
    rx_nwchem_version = regex(r'^\s*Northwest Computational Chemistry Package \(NWChem\)\s+([-0-9a-zA-Z_.]+)')
    rx_e_scf = regex(r'^\s*Total (SCF|DFT) energy *= *([-+eE0-9.]+)')
    rx_e_ccsd = regex(r'^\s*CCSD total energy */ *hartree *= *([-+eE0-9.]+)')
    rx_e_ccsd_t = regex(r'^\s*CCSD\(T\) total energy */ *hartree *= *([-+eE0-9.]+)')
    rx_e_nucl = regex(r'^\s*Nuclear repulsion energy *= *([-+eE0-9.]+)')

    txtfile = text_input(filename)

    self.clear()
    rslt = self

    # This will also serve as an initial screening of the file
    try:
      L = txtfile.seek_text(rx_nwchem_version.rx)
    except:
      raise RuntimeError, \
        "Cannot determine nwchem version in `%s'; perhaps it is not an nwchem output file" % filename
    rslt['info_code_version'] = (rx_nwchem_version % L).group(1)
    rslt['info_mtime'] = time.localtime(os.path.getmtime(filename))

    search_patterns = [
      (regex(r'^\s*Total SCF energy *= *([-+eE0-9.]+)'),                   'E_SCF', float),
      (regex(r'^\s*Total DFT energy *= *([-+eE0-9.]+)'),                   'E_DFT', float),
      (regex(r'^\s*CCSD total energy */ *hartree *= *([-+eE0-9.]+)'),      'E_CCSD', float),
      (regex(r'^\s*CCSD\(T\) total energy */ *hartree *= *([-+eE0-9.]+)'), 'E_CCSD_T', float),
      (regex(r'^\s*Nuclear repulsion energy *= *([-+eE0-9.]+)'),           'E_nuclear', float),
      (regex(r'^\s*Total MP2 energy\s+([-+eE0-9.]+)'),                     'E_MP2', float), # old MP2 module
    ]

    for L in txtfile:
      flds = L.split()
      if len(flds) == 0:
        continue
      else:
        for (pat, act, arg1) in search_patterns:
          if pat % L:
            if isinstance(act, str):
              rslt[act] = arg1(pat[1])
              break

    return rslt

