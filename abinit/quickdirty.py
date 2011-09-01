# $Id: quickdirty.py,v 1.2 2011-09-01 20:11:33 wirawan Exp $
#
# pyqmc.abinit.quickdirty
#
# Created: 20100210
# Wirawan Purwanto
#
# Quick and dirty tools for ABINIT

import sys
import os
import os.path

from wpylib.file.file_utils import open_input_file

def abinit_get_bands(file):
  """Reads in ABINIT band information (eigenvalue, occupancy) from an output
  file.
  We assuming that there is only one ABINIT calculation per output file,
  with occupancy data printed right after the eigenvalues."""
  F = open_input_file(file)
  up = []
  dn = []
  for ll in F:
    L = ll.strip()
    flds = L.split()
    if L.startswith("Eigenvalues (hartree) for nkpt"):
      nkpts = int(flds[4])
      if L.find("SPIN DOWN") >= 0:
        band_sector = dn
      else:
        band_sector = up

      for ikpt in xrange(nkpts):
        flds = F.next().replace("="," ").split()
        if flds[0] != "kpt#": raise AssertionError, "Expecting k-point info here"
        nbands = int(flds[3].rstrip(","))
        kpt = tuple([ float(flds[i]) for i in [7,8,9] ])
        iband = 0
        egnvals = []
        while iband < nbands:
          flds = F.next().split()
          egnvals += [ float(f) for f in flds ]
          iband += len(flds)
        if not F.next().strip().startswith("occupation numbers for kpt#"):
          raise AssertionError, "Expecting occupation number info here"
        occ = []
        iband = 0
        while iband < nbands:
          flds = F.next().split()
          occ += [ float(f) for f in flds ]
          iband += len(flds)
        band_sector.append({
            'nbands' : nbands,
            'kpt' : kpt,
            'egnvals' : egnvals,
            'occ' : occ,
          })
  F.close()
  return (up, dn)



