# -*- python -*-
#
# pyqmc.utils.cost module
#
# Wirawan Purwanto
# Created: 20130514
#
# Based on $PWQMC77/scripts/cost.py, with CVS revision info:
#    v 1.10 2013-03-27 15:18:01 wirawan Exp
#
#

"""
pyqmc.utils.cost

AFQMC cost estimator and analyzer.

"""

import math
import os
import sys

from math import sqrt
from wpylib.iofmt.text_output import text_output


class cost_base(object):
  """Crude base class for cost estimator."""
  def def_value(self, name, value):
    """Sets a default value, or return the existing value."""
    if value == None:
      value = getattr(self, name)
      return value
    else:
      setattr(self, name, value)
      return value
  def setvals(self, **args):
    """A quick way to set multiple attributes."""
    for (k,v) in args.iteritems():
      setattr(self, k, v)



class qmc_cost_estimator(cost_base):
  """QMC cost calculator, given the number of blocks etc."""
  size_complex = 16
  size_real = 8

  @property
  def numprocs(self): # old keyword
    return self.num_tasks

  """The 'get_XXX' mechanism below allows us to define reasonable
  defaults if the XXX property is not set."""

  @property
  def get_nblkstep(self):
    if not hasattr(self, "nblkstep"):
      nblkstep = int(round(self.betablk / self.deltau))
    else:
      nblkstep = self.nblkstep
    return nblkstep

  @property
  def get_nwlkmin(self):
    return getattr(self, "nwlkmin", max(self.nwlk//2, 1))

  @property
  def get_nwlkmax(self):
    return getattr(self, "nwlkmax", 2*self.nwlk)

  @property
  def get_nwlkmax_proc(self):
    nwlkmax_proc = getattr(self, "nwlkmax_proc",
                           (self.get_nwlkmax + self.num_tasks - 1) // self.num_tasks)
    return nwlkmax_proc

  @property
  def get_hsop_dim(self):
    if getattr(self, "general_H2", False):
      return self.nbasis * (self.nbasis + 1) // 2
    else:
      return self.nbasis**2

  def compute_time_cost(self, Print=False):
    """
    Get total computing time cost.

    Required input:
    - deltau
    -
    """
    nblkstep = self.get_nblkstep
    #if not hasattr(self, "nblkstep"):
    #  nblkstep = int(round(self.betablk / self.deltau))
    #else:
    #  nblkstep = self.nblkstep
    self.neqstep = self.neq * nblkstep
    self.ngrthstep = self.ngrth * nblkstep
    self.nmeastep = self.nblk * nblkstep
    self.beta_eq = self.neqstep * self.deltau
    self.beta_grth = self.ngrthstep * self.deltau
    self.beta_meas = self.nmeastep * self.deltau

    nwlkmax = getattr(self, "nwlkmax", 2*self.nwlk)
    #nwlkmax_proc = getattr(self, "nwlkmax_proc",
    #                       (nwlkmax + self.num_tasks - 1) // self.num_tasks)
    nwlkmax_proc = self.get_nwlkmax_proc

    self.nblkTotal = self.neq + self.ngrth + self.nblk
    self.nwlkstepmax_proc_tot = self.nblkTotal * nwlkmax_proc # max num of wlk-steps per MPI task

    # Wall times begin with "t"
    # FIXME: for now we ignore popctl and modgs; later these have to be accounted for.
    self.teq = self.neqstep * nwlkmax_proc * self.timing['step']
    self.tgrth = self.ngrthstep * nwlkmax_proc * self.timing['step']
    self.tmeas1 = nblkstep * nwlkmax_proc * self.timing['step']
    self.tmeas = self.nmeastep * nwlkmax_proc * self.timing['step']

    # Total computer-hour times begin with "T"
    # computer-hour metric:
    self.Teq = self.teq * self.num_tasks * self.num_threads
    self.Tgrth = self.tgrth * self.num_tasks * self.num_threads
    self.Tmeas1 = self.tmeas1 * self.num_tasks * self.num_threads
    self.Tmeas = self.tmeas * self.num_tasks * self.num_threads

    if Print: self.printout_timing()

  def compute_mem_cost(self, Print=False):
    """
    Get total computing memory cost.

    Required input:
    - nbasis
    - nflds
    - nwlkmax
    - nptot, nup, ndn
    - npsitdet
    """
    raise NotImplementedError, "Under construction"
    nwlkmax_proc = self.get_nwlkmax_proc
    (M, Nptot, F, D) = (self.nbasis, self.nptot, self.nflds, self.npsitdet)
    (DPC, DP) = (self.size_complex, self.size_real)
    """ UNDER CONSTRUCTION
    self.wlk_size = Nptot * M * DPC
    self.mem_wlk = self.wlk_size * nwlkmax_proc
    self.mem_Lvec = self.get_hsop_dim * F * DP   # so far it is double precision
    self.mem_
    """


    if Print: raise NotImplementedError

  def printout_mem(self, out=None):
    """
    Prints out a report for memory estimate.
    """


  def printout_timing(self, out=None):
    """
    Prints out a report for timing estimate.
    """
    if out == None: out = sys.stdout
    O = text_output(out, flush=True)
    #info = self.info
    nblkstep = self.get_nblkstep
    h = 60 * 60  # secs per hour
    Mh = h * 1.0e6  # million secs per hour
    O("deltau = %.14g\n" % self.deltau)
    O("\n")
    O("beta per blk = %10.3f  ( %9d steps = %10.4f M avg wlk-steps )\n" % (nblkstep*self.deltau, nblkstep, nblkstep*self.nwlk*1e-6))
    O("beta_eq      = %10.3f  ( %9d steps = %10.4f M avg wlk-steps )\n" % (self.beta_eq, self.neqstep, self.neqstep*self.nwlk*1e-6))
    O("beta_grth    = %10.3f  ( %9d steps = %10.4f M avg wlk-steps )\n" % (self.beta_grth, self.ngrthstep, self.ngrthstep*self.nwlk*1e-6))
    O("beta_meas    = %10.3f  ( %9d steps = %10.4f M avg wlk-steps = %.6f M avg wlk-beta )\n"
      % (self.beta_meas, self.nmeastep,
         self.nmeastep*self.nwlk*1e-6,
         self.nmeastep*self.deltau*self.nwlk*1e-6))
    O("beta_total   = %10.3f  ( %9d steps = %10.4f M avg wlk-steps )\n" % (
        self.beta_eq + self.beta_grth + self.beta_meas,
        self.neqstep + self.ngrthstep + self.nmeastep,
        (self.neqstep + self.ngrthstep + self.nmeastep) * self.nwlk * 1e-6,
    ))
    O("\n")
    O("num_tasks    = %d\n" % (self.num_tasks))
    O("num_threads  = %d\n" % (self.num_threads))
    O("num_cores    = %d total\n" % (self.num_tasks * self.num_threads))
    O("\n")
    O("nwlkmin      = %d\n" % (self.get_nwlkmin))
    O("nwlk         = %d\n" % (self.nwlk))
    O("nwlkmax      = %d\n" % (self.get_nwlkmax))
    O("nwlkmax_proc = %d\n" % (self.get_nwlkmax_proc))
    O("\n")
    O("wallclock stats:\n")
    O("teq          = %10.3f h\n" % (self.teq / h))
    O("tgrth        = %10.3f h\n" % (self.tgrth / h))
    O("tmeas1       = %10.3f h (%.2s secs) per measurement blk\n" % (self.tmeas1 / h, self.tmeas1))
    O("tmeas        = %10.3f h\n" % (self.tmeas  / h))
    O("ttotal       = %10.3f h\n" % ((self.teq + self.tgrth + self.tmeas) / h))
    O("\n")
    O("overall computer resource stats:\n")
    O("Teq          = %10.6f M core*h\n" % (self.Teq / Mh))
    O("Tgrth        = %10.6f M core*h\n" % (self.Tgrth / Mh))
    O("Tmeas1       = %10.6f M core*h\n" % (self.Tmeas1 / Mh))
    O("Tmeas        = %10.6f M core*h\n" % (self.Tmeas / Mh))
    O("Ttotal       = %10.6f M core*h\n" % ((self.Teq + self.Tgrth + self.Tmeas) / Mh))



class real_qmc_cost(object):
  """A class to parse in GAFQMC or PWQMC INFO file and analyze the
  real computational cost of that calculation.
  """
  def read(self, filename, clock_shift=0):
    """
    Optional argument clock_shift defines the clock shift for
    the start and end LOCAL times of the computer where the calculation
    took place.
    For example, if your local time is EST, and the computer where you did the
    calculation has CST clock, then clock_shift is -1.
    We will adjust the start and end times to match your local time.
    """
    from pyqmc.results.gafqmc_info import gafqmc_info
    from pyqmc.utils.gafqmc_quick_dirty import is_qmc_output as is_gafqmc_output
    from pyqmc.results.pwqmc_info import pwqmc_info
    from wpylib.datetime import shift_time
    if is_gafqmc_output(filename):
      self.info = gafqmc_info(filename)
    else:
      self.info = pwqmc_info(filename)
    info = self.info
    if clock_shift != 0:
      dt = -clock_shift * 3600.0
      if 'start_time' in info:
        info['start_time'] = shift_time(info['start_time'], dt)
      if 'end_time' in info:
        info['end_time'] = shift_time(info['end_time'], dt)
    if not "num_threads" in info:
      info['num_threads'] = 1
    if not "num_tasks" in info:
      info['num_tasks'] = 1


  def analyze(self):
    info = self.info
    self.walltime = time_diff(info["end_time" if "end_time" in info else "info_mtime"], info.start_time) # in secs
    self.cputime = self.walltime / 3600.0 * info.num_tasks * info.num_threads # in core-hours
    self.nblk_actual = len(info.meas_energy)
    self.nsteps_all = (info.neq + info.ngrth + self.nblk_actual) * info.nblkstep  # total QMC population-step count
    self.nwlkmax_proc = (info.nwlkmax + info.num_tasks - 1) // info.num_tasks
    self.nwlkmin_proc = (info.nwlkmin) // info.num_tasks
    self.nwlkavg_proc = float(info.nwlk) / info.num_tasks

    # nwsteps = number of walker-steps (total during QMC run), per MPI task
    self.nwsteps_min = self.nsteps_all * self.nwlkmin_proc   # est. by nwlkmin_proc
    self.nwsteps_avg = self.nsteps_all * self.nwlkavg_proc   # est. by nwlkavg_proc
    self.tsteps_max = self.walltime / self.nwsteps_min
    self.tsteps_avg = self.walltime / self.nwsteps_avg


  def printout(self, out=None):
    from pyqmc.results.gafqmc_info import gafqmc_info
    from pyqmc.results.pwqmc_info import pwqmc_info
    from wpylib.iofmt.text_output import text_output
    if out == None: out = sys.stdout
    O = text_output(out, flush=True)
    info = self.info
    if isinstance(info, gafqmc_info):
      O("Calculation: GAFQMC\n")
    elif isinstance(info, gafqmc_info):
      O("Calculation: PWQMC\n")

    beta_meas = self.nblk_actual * self.info.betablk
    O("walltime          = %g secs  (%g hours)\n" % (self.walltime, self.walltime/3600.0))
    O("cpu core count    = %d\n" % (info.num_tasks * info.num_threads))
    O("cputime          >= %.6f M core*h\n" % (self.cputime * 1e-6))
    O("num meas blocks   = %d    ( beta = %.3f  <-->  %.4f M avg wlk-steps = %.6f M avg wlk-beta )\n"
      % (self.nblk_actual, beta_meas,
         self.nblk_actual * info.nblkstep * info.nwlk * 1e-6,
         beta_meas * info.nwlk * 1e-6,
        ))
    O("num all blocks    = %d\n" % (info.neq + info.ngrth + self.nblk_actual))
    O("num all steps     = %d\n" % (self.nsteps_all))
    O("nwlkmax_proc      = %d\n" % (self.nwlkmax_proc))
    O("\n")
    O("AVERAGE CASE\n")
    O("nwlkavg_proc      = %.10g\n" % (self.nwlkavg_proc))
    O("avg tot wlk steps = %.10g (per MPI task)\n" % (self.nwsteps_avg))
    O("avg step time     = %.10g\n" % (self.tsteps_avg))
    O("\n")
    O("WORST CASE\n")
    O("nwlkmin_proc      = %.10g\n" % (self.nwlkmin_proc))
    O("min tot wlk steps = %.10g (per MPI task)\n" % (self.nwsteps_min))
    O("WORST step time   = %.10g\n" % (self.tsteps_max))


def time_diff(time1, time2):
  from time import mktime
  return mktime(time1) - mktime(time2)
