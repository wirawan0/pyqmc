#!/usr/bin/env python
#
# Wirawan Purwanto
# Created: 20140505
#

"""
pyqmc.bin.convert_BSE_basis

A utility to convert BSE basis text (in GAMESS format) to
an internal format for pyqmc use.

The input must be in GAMESS-US format.
"""

import sys

from pyqmc.basis.gaussian import convert_BSE_basis

def main(argv):
  """Main program
  """
  inpfile = argv[0]
  try:
    outfile = argv[1]
  except IndexError, KeyError:
    outfile = sys.stdout
  convert_BSE_basis(inpfile, outfile)


if __name__ == "__main__":
  main(sys.argv[1:])

