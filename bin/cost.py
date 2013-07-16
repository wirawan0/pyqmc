#!/usr/bin/env python
#
# cost.py
#
# Wirawan Purwanto
# Created: 20130220
#

"""
pyqmc.bin.cost

Utility to analyze or estimate AFQMC calculation cost.
"""

import os
import sys
from pyqmc.utils.cost import real_qmc_cost

def Help():
  help_str = """\
Available commands:

  actual  <INFO_file>  [clock_shift]

    Computes the actual cost of an AFQMC calculation, especially,
    estimating the amount of time required to propagate a walker.
"""


def analyze_real_calc_cost(filename, **_opts_):
  """One-shot function for cost analysis of a real AFQMC calculation.
  """
  opts = _opts_
  Cost = real_qmc_cost()
  Cost.read(filename, **opts)
  Cost.analyze()
  Cost.printout()
  return Cost



def main():
  if len(sys.argv) <= 1:
    print >>sys.stderr, "Nothing to do?"
    return 1
  elif sys.argv[1] in ('actual', '--actual', '--real'):
    filename = sys.argv[2]
    if len(sys.argv) > 3:
      clock_shift = float(sys.argv[3])
    else:
      clock_shift = 0
    analyze_real_calc_cost(filename, clock_shift=clock_shift)
  else:
    print >>sys.stderr, "Invalid command:", sys.argv[1]
    return 1


if __name__ == "__main__":
  sys.exit(main())

