# $Id: pwqmc_info.py,v 1.2 2010-08-13 01:55:18 wirawan Exp $
#
# pwqmc_info.py
# Tools to parse PWQMC-77 INFO file
#
# Wirawan Purwanto
# Created: 20100409
#

# Python standard modules
import math
import os
import sys
import time

import numpy

class pwqmc_info(dict):
  '''Structure to represent the metadata contained in INFO file.

  Available information:
  * info_file
  * Evar_noconst  Evar  H0
  * deltau  betablk
  * kpt
  * vol
  '''
  def __init__(self, src=None):
    if isinstance(src, dict):
      self.clear()
      self.update(src)
    if isinstance(src, str):
      self.parse_INFO(src)
      self.filename = src
    else:
      pass # self._data = {}
  def __getattr__(self, key):
    try:
      return self[key]
    except:
      return dict.__getattribute__(self, key)
  # __setattr__ is not set. INFO-related attributes are supposed to be
  # read-only.
  def parse_INFO(self, INFO):
    '''Gets all the necessary info (calculation parameters) from the INFO file.
    This is a very old routine.
    We use this as temporary starting point.'''
    info_file = open(INFO, "r")
    self.clear()
    rslt = self
    rslt['info_file'] = INFO
    for L in info_file:
      Ls = L.strip()
      flds = Ls.split()
      if len(flds) == 0:
        continue
      elif flds[0] == "Subtotal":
        rslt["Evar_noconst"] = float(flds[2])
      elif flds[0] == "Variational" and flds[1] == "energy":
        rslt["Evar"] = float(flds[3])
        rslt["H0"] = rslt["Evar"] - rslt["Evar_noconst"]
      elif flds[0] == "deltau,":
        rslt["deltau"] = float(flds[3])
      elif flds[0] == "beta=":
        rslt["betablk"] = float(flds[1])
      elif Ls.startswith("Using reduced k-pts:"):
        kx = float(flds[3])
        ky = float(flds[4])
        if len(flds) > 5:
          kz = float(flds[5])
        else:
          kz = float(info_file.next().split()[0])
        rslt["kpt"] = (kx,ky,kz)
      elif Ls.startswith("Volume of the unit cell ="):
        rslt["vol"] = float(flds[6]) # in bohr**3
      elif Ls.startswith("Input Etrial="):
        rslt["Etrial_noconst"] = float(flds[2]) # no H0 yet
        #print Ls
      elif Ls.startswith("New etrial to be used in El_bound:"):
        rslt["Etrial_noconst"] = float(flds[7]) # no H0 yet
        #print Ls
      elif Ls.startswith("read in new anorm + etrial:"):
        rslt["anorm"] = float(flds[6])
        if len(flds) > 7:
          rslt["Etrial_noconst"] = float(flds[7])
        else:
          rslt["Etrial_noconst"] = float(info_file.next().split()[0])
        #print Ls
      elif Ls.startswith("itv_Em="):
        rslt["itv_Em"] = int(flds[1])
      elif Ls.startswith("itv_pc="):
        rslt["itv_pc"] = int(flds[1])
      elif Ls.startswith("itv_pc_eq="):
        rslt["itv_pc_eq"] = int(flds[1])
      elif Ls.startswith("nblkstep="):
        rslt["nblkstep"] = int(flds[1])
      elif Ls.startswith("nwlk="):
        rslt["nwlk"] = int(flds[1])
      elif Ls.startswith("nwlkmax="):
        rslt["nwlkmax"] = int(flds[1])
      elif Ls.startswith("nwlkmin="):
        rslt["nwlkmin"] = int(flds[1])
    info_file.close()
    rslt.setdefault("nwlkmax", rslt.nwlk * 2)
    rslt.setdefault("nwlkmin", max(rslt.nwlk / 2, 1))
    return rslt
  def __str__(self):
    return "<" +self.__module__ + "." + self.__class__.__name__ + \
      " object (" + dict.__str__(self) + ")>"



# VERY VERY OLD API: used as starting point

def parse_INFO(INFO):
  '''Gets all the necessary info (calculation parameters) from the INFO file.'''
  return pwqmc_info(str(INFO))


def kptstr(kpt):
  """Prints k-point string in a standardized way.
  kpt is either a string (in which case it is passed verbatimly), or
  a 3-vector (list, tuple, numpy array, whatever indexable with integers
  0, 1, 2)."""
  if isinstance(kpt, str):
    # FIXME: must check assumption
    # But this is also nice for "*" kind of wildcard.
    # Assume already +3425+3425+3425 format:
    return kpt
  else:
    return "%+05.0f%+05.0f%+05.0f" \
           % (float(kpt[0])*10000, float(kpt[1])*10000, float(kpt[2])*10000)



