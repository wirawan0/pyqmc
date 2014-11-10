"""Test codes for gafqmc_info.
"""

import os
import os.path
import numpy

from pprint import pprint


def SPECIMEN_FILE(fn):
  return os.path.join("SPECIMENS", fn)

def test_gafqmc_info_BugStop_in_meas():
  """[20140905] Tests the ability of result parser to recover from BugStop-ed
  run output.
  """
  from pyqmc.results.gafqmc_info import gafqmc_info
  File = SPECIMEN_FILE("INFO.BugStop-in-meas")
  Info = gafqmc_info(File)
  print "INFO output:"
  pprint(Info)
  if "BUGSTOP" in Info:
    print "test_gafqmc_info_BugStop_in_meas: OK"
  else:
    print "test_gafqmc_info_BugStop_in_meas: FAIL"

if __name__ == "__main__":
  test_gafqmc_info_BugStop_in_meas()

