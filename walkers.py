# $Id: walkers.py,v 1.1 2009-01-09 22:06:11 wirawan Exp $
#
# walkers.py
# Created: 20080930
# Wirawan Purwanto
#
# Handling of AFQMC walker file in Python.

import sys
import os
import subprocess  # for piping
from pyqmc.multidet import Det, MultiDet, mdet_ovlp
from pyqmc.matrices import read_det_matrix

#class AFQMCPopMetadata(object):
  #def __init__(self):
  #  pass

class AFQMCPop(MultiDet):
  # Static vars:
  meta_list = {
    "nbasis" : int,
    "nptot" : int,
    #"lran" : int(4), # disable this for now
    "nblktot" : int,
    "uptot" : float,
    "downtot" : float,
    "srun" : float,
    "s2run" : float,
    "fmt_version" : float,
    "nh" : int,
    "anorm" : float,
    "etrial" : float,
    "istpacc" : int,
    "timeshift" : float,
    "nwlkr_proc" : int,
    "iflg_chkpt_impfn" : int,
    # for now these are kept as string:
    "day" : str,
    "time" : str,
  }
  analyze_walker_script = os.environ["GAFQMC"] + "/scripts/analyze-walkers"
  # Members:
  # dets = list of walkers
  # other members are listed above in meta_list
  def __init__(self, fname = None, PsiT = None, restricted = None,
               nup = None, ndn = None, verbose = None):
    if fname != None:
      self.nup = nup
      self.ndn = ndn
      self.read(fname, PsiT, restricted, verbose)

  def check_nparts(self):
    if (self.nup == None or self.ndn == None):
      raise ValueError, \
        "Number of particles (nup, ndn) has not been initialized properly."

  #def analyze_walker_script(self):
  #  return os.environ["GAFQMC"] + "/scripts/analyze-walkers"

  def read(self, fname, PsiT = None, restricted = None, verbose = None):
    self.read_gafqmc_pop(fname, PsiT, restricted, verbose)

  def read_gafqmc_pop(self, fname, PsiT = None,
                      restricted = None, verbose = None): #{
    '''Reads in GAFQMC population into a MultiDet-style wave function.
    '''
    self.check_nparts()
    nup = self.nup
    ndn = self.ndn

    meta_list = AFQMCPop.meta_list
    cmnd = [ self.analyze_walker_script, fname, "--dumpwlkrs" ]

    # Always start afresh: delete all previously set data
    for k in AFQMCPop.meta_list.keys():
      if (hasattr(self, k)): delattr(self, k)

    self.dets = []
    #self.phase = []
    #self.El = []
    #self.impfn = []
    if (verbose): print "Executing: ", " ".join(cmnd)
    px = subprocess.Popen(cmnd, shell=False, bufsize=-1,
                          stdout=subprocess.PIPE)
    pipe = px.stdout
    s = pipe.readline()
    while s != "":
      s = s.strip()
      if s.startswith("="):
        # Special exception devised for "========" kind of separator.
        pass
      elif "=" in s:
        #print "px: ", s
        [ kwd, val ] = s.split("=", 1)
        kwd = kwd.strip()
        val = val.strip()
        if kwd in meta_list:
          # FIXME: Right now we only support scalar attributes.
          #print kwd, meta_list[kwd],
          setattr(self, kwd, map(meta_list[kwd], [val])[0])
          if (verbose):
            print "set: ", kwd, "->", val, " == ", getattr(self, kwd)
      elif s.startswith("Wlk#"):
        # next subloop to read in the walkers
        for i in xrange(0, self.nwlkr_proc):
          s = pipe.readline()
          fld = s.split()
          #d = Det()
          d = read_det_matrix(pipe, "walker#" + fld[0],
                              self.nbasis, self.nup, self.ndn,
                              cplx=True, restricted=restricted,
                              verbose=verbose)
          d.weight = float(fld[1])
          d.phase = float(fld[2])
          d.El = float(fld[3]) + 1j*float(fld[4])
          self.dets.append(d)
      elif s.startswith("impfn data:"): # This is <psiT|phi>
        if (verbose): print "Reading impfn data"
        for d in self.dets:
          s = pipe.readline().strip(" ()\n\r\t")
          #print s
          [s1, s2] = s.split(",")
          d.impfn = complex(float(s1), float(s2))
          d.ampl = d.weight / d.impfn
      #print ">", s
      s = pipe.readline()

    pipe.close()
    px.wait()

    if not hasattr(self.dets[0], "impfn"):
      if PsiT != None:
        rdet = MultiDet()
        rdet.dets = [ 0 ] # dummy init
        for d in self.dets:
          rdet.dets[0] = d
          d.ampl = 1.0
          d.impfn = mdet_ovlp(PsiT, rdet)
          d.ampl = d.weight / d.impfn
      else:
        sys.stderr.write("Warning: cannot deduce the amplitude of " + \
                         "GAFQMC population in file `" + fname + \
                         "' because PsiT argument is not given.\n")
  #}read_gafqmc_pop
