# $Id: constants.py,v 1.1 2010-09-07 15:09:35 wirawan Exp $
#
# pyqmc.libs.constants module
# Created: 20100902
# Wirawan Purwanto
#
"""
Library of constants from various computational packages.
At present it contain constants from:
* nwchem 5.1.1
* MPQC 2.3.1
"""


class nwchem_5_1_1:
  """Constants from NWCHEM 5.1.1 .

  References: geom.F (where the printout of the geometries is done)
  """

  # The default scaling is 1.88972598858 a.u. per angstrom
  # defined in geom_input.F
  ang2bohr = 1.88972598858
  #bohr_in_ang = 0.52917715 -- in geom_hnd.F, perhaps not used for
  # this purpose
  bohr2ang = 1.0 / ang_in_bohr

class mpqc_2_3_1:
  """Constants in MPQC 2.3.1 .

  References: mpqc-2.3.1/src/lib/util/misc/units.cc
  """

  # Old from CRC Handbook 77th Ed. Tables 1&2 (CODATA 1986)
  # But at least, it is consistent with nwchem 5.1.1 constant:
  bohr2ang = 5.29177249e-01
  ang2bohr = 1.0 / bohr2ang
  hartree2eV = 27.2113834

