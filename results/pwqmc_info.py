# $Id: pwqmc_info.py,v 1.4 2010-11-11 18:03:46 wirawan Exp $
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
import re
import sys
import time

import numpy

from wpylib.iofmt.text_input import text_input
from wpylib.db.result_base import result_base
from wpylib.sugar import ifelse

class pwqmc_info(result_base):
  '''Structure to represent the metadata contained in INFO file.

  Available information:
  * info_file
  * trial_wfn_file
  * nelec_up  nelec_dn  udet
  * nup  ndn   (deprecated name of nelec_up and nelec_dn)
  * nbasis  LL[0:3]
  * Evar_noconst  Evar  H0
  * Etrial_noconst  Etrial
  * deltau  betablk  nblkstep
  * nwlk nwlkmax nwlkmin
  * itv_Em itv_pc itv_pc_eq
  * kpt
  * vol

  Additional runtime info will also be parsed in:

  '''

  meas_dtype = numpy.dtype([('beta',float), ('overlap',float),
                            ('Etotal',float), ('Eproj',float)])
  # WARNING: Jellium is weird (FIXME, TODO)--Eproj and Etotal are
  # reversed, because only the correlation energy is what people care
  # here.
  meas_jell_dtype = numpy.dtype([('beta',float), ('overlap',float),
                                 ('Eproj',float), ('Etotal',float)])

  def parse_INFO(self, INFO):
    '''Gets all the necessary info (calculation parameters) from the INFO file.
    This is a very old routine.
    We use this as temporary starting point.'''
    # FIXME: comment_char is temporarily set to ASCII 0, which should be
    # an invalid character in this output file.
    info_file = text_input(INFO, comment_char='\0', skip_blank_lines=False)
    self.clear()
    rslt = self
    rslt['info_file'] = INFO
    rslt['info_mtime'] = time.localtime(os.stat(INFO).st_mtime)
    for L in info_file:
      Ls = L.strip()
      ls = Ls.lower()
      flds = Ls.split()
      if len(flds) == 0:
        continue
      elif Ls.startswith("# of particles:"):
        u = int(flds[3])
        d = int(flds[4])
        if u < d:
          sys.stderr.write(
            "pwqmc_info.parse_INFO:Warning: nup < ndn in info file `%s'; autofixing this mistake!\n" % (INFO)
          )
          t = u; u = d; d = t
        rslt['nup'] = u
        rslt['ndn'] = d
        rslt['nelec_up'] = u
        rslt['nelec_dn'] = d
      elif ls.startswith("majority and minority det are coupled"):
        rslt['udet'] = False
      elif ls.startswith("majority and minority det are independent"):
        rslt['udet'] = True
      elif flds[0] == "Nbasis":
        rslt['nbasis'] = int(flds[2])
      elif ls.startswith("input fft dimension ll ="):
        rslt['LL'] = (int(flds[5]), int(flds[6]), int(flds[7]))
      elif ls.startswith("trial wf from input: "):
        rslt['trial_wfn_file'] = Ls[20:].strip()
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

      # ---crystal and atom parameters, etc---
      # Wherever possible, we will use abinit-style keywords
      # to make this consistent with abinit.
      elif Ls.startswith("Volume of the unit cell ="):
        rslt["vol"] = float(flds[6]) # in bohr**3
      elif Ls.startswith("nspec_in ="):
        rslt["ntypat"] = int(flds[2]) # number of atomic species
        rslt["has_pseudo"] = (rslt["ntypat"] > 0)

      # ---runtime info below---
      elif Ls.startswith("OpenMP support enabled using"):
        rslt["num_threads"] = int(flds[4])
      elif Ls.startswith("parallel version using"):
        rslt["num_tasks"] = int(flds[3])
      elif Ls.startswith("hostname:"):
        rslt["run_host"] = flds[1]
      elif Ls.startswith("Program started on"):
        rslt["start_time"] = \
          time.strptime(flds[3] + " " + flds[5][:6], "%Y%m%d %H%M%S")

      # measurement and other complex data capture
      elif Ls.startswith("Measurement phase..."):
        self.parse_measurement0(info_file, rslt)

    info_file.close()
    rslt.setdefault("nwlkmax", rslt.nwlk * 2)
    rslt.setdefault("nwlkmin", max(rslt.nwlk / 2, 1))
    rslt["run_mpi"] = ("num_tasks" in rslt)
    rslt["run_openmp"] = ("num_threads" in rslt)
    return rslt
  parse_file_ = parse_INFO


  def parse_measurement0(self, info_file, rslt):
    """Internal routine to parse only the measurement results of the
    file.
    info_file is an open file-like object.
    The last line read must have been 'Measurement phase...'

    TODO:
    - add stand-alone parse_measurement routine?
    """
    from pyqmc import PyqmcParseError

    L1 = info_file.next().strip()
    L2 = info_file.next().strip()
    if not L1.startswith("psi*psi_T") or not L2.startswith("-----"):
      raise PyqmcParseError, \
        "Error detecting the beginning of measurement phase output stream."

    for_D2E = lambda s : s.replace("D","E").replace("d","e")
    EOS = re.compile(r"^\s*\*+\s*$")  # end-of-stream marker
    RS = re.compile(r"^\s*-+\s*$")  # record separator

    meas = []
    for L in info_file:
      Ls = for_D2E(L.strip())
      flds = Ls.split()
      if EOS.search(Ls):
        break   # end-of-stream detected
      elif len(flds) == 3:
        # special case to handle wrapped Fortran output
        Ls2 = for_D2E(info_file.next().strip())
        flds2 = Ls2.split()
        if len(flds2) == 0:
          raise PyqmcParseError, \
            "Invalid format in PWQMC measurement text (INFO)"
        flds.append(flds2[0])
      elif len(flds) < 4:
        raise PyqmcParseError, \
          "Invalid format in PWQMC measurement text (INFO)"
      try:
        rec = tuple(map(float, flds[:4]))
      except:
        raise PyqmcParseError, \
          "Error parsing PWQMC measurement text (INFO)"
      meas.append(rec)

      Ls = info_file.next().strip()
      if not RS.search(Ls):
        raise PyqmcParseError, \
          "Invalid record separator in PWQMC measurement text (INFO)"

    dtype = ifelse(rslt["has_pseudo"], self.meas_dtype, self.meas_jell_dtype)
    rslt["meas_energy"] = numpy.array(meas, dtype=dtype)



# VERY VERY OLD API: used as starting point

def parse_INFO(INFO):
  '''Gets all the necessary info (calculation parameters) from the INFO file.'''
  return pwqmc_info(str(INFO))


def kptstr(kpt):
  """Prints k-point string in a standardized way.
  kpt is either a string (in which case it is passed verbatimly), or
  a 3-vector (list, tuple, numpy array, whatever indexable with integers
  0, 1, 2)."""
  if isinstance(kpt, basestring):
    # FIXME: must check assumption
    # But this is also nice for "*" kind of wildcard.
    # Assume already +3425+3425+3425 format:
    return kpt
  else:
    return "%+05.0f%+05.0f%+05.0f" \
           % (float(kpt[0])*10000, float(kpt[1])*10000, float(kpt[2])*10000)



