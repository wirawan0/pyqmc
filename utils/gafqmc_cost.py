# -*- python -*-
#
# pyqmc.utils.gafqmc_cost module
#
# Wirawan Purwanto
# Created: 20130923
#
#

"""
pyqmc.utils.gafqmc_cost

Cost estimator and analyzer for GAFQMC code.

"""

import numpy

from pyqmc.utils import cost
from pyqmc.utils import linalg_cost


class gafqmc_sparse1_cost_estimator(cost.qmc_cost_estimator):
  """Cost estimator specifically for GAFQMC calculation.

  Names of precomputed objects:
  * Vijkl = four-indexed two-body operator
  * Vss = product of two trial wfn orbitals (same spin) with Vijkl
  * Vxs = product of two trial wfn orbitals (opposite spins) with Vijkl
  * Ls  = product of one trial wfn orbital (for all spins) with L
          one-body operator.
  """

  # (SPARSE) PRECOMPUTED MATRICES

  # Default sparse matrix density (ballpark estimate)
  # These values MAY NOT BE CORRECT for your particular calculation!
  # These are ok for GTO-basis calculations without frozen core:
  dens_Vss = 0.5
  dens_Vxs = 0.5
  dens_Ls = 0.7

  tpref_gf1_ovlp = 9.74193418e-09
  tpref_gf1_ovlpinv = 2.25385979e-09
  tpref_FB = 1.83499926e-08
  tpref_Elocal = 1.09604036e-08


  # Turns out, empirically Vss and Vxs densities are the same

  def __init__(self):
    self.linalg = linalg_cost.linalg_cost_ops()

  def compute_mem_cost(self, Print=False):
    """
    Estimate a calculation's MEMORY cost based on the given input sizes.
    For original sparse method due to Wissam.

    Required input:
    - nbasis
    - nflds
    - nwlkmax
    - nptot, nup, ndn
    - npsitdet

    Objects:

    Output:
    - mem_Vss =
    - mem_Vxs =
    - mem_rhos =
    """
    nwlkmax_proc = self.get_nwlkmax_proc
    (M, Nptot, Nu, Nd, F, D) = self.params_wlkr(0)
    (dpc, dp, it) = self.params_sys(0)

    self.wlk_size = Nptot * M * dpc
    self.mem_wlk = self.wlk_size * nwlkmax_proc
    self.mem_Lvec = self.get_hsop_dim * F * dp   # so far it is double precision
    # number of elements in the sparse objects, per determinant
    self.count_Vuu_det = self.dens_Vss * M**2 * Nu**2
    self.count_Vdd_det = self.dens_Vss * M**2 * Nd**2
    self.count_Vud_det = self.dens_Vxs * M**2 * Nu*Nd
    self.count_Ls_det = self.dens_Ls * F * M * (Nu+Nd)
    # number of elements in the sparse objects, ALL determinants
    self.count_Vuu = D * self.count_Vuu_det
    self.count_Vdd = D * self.count_Vdd_det
    self.count_Vud = D * self.count_Vud_det
    self.count_Ls = self.dens_Ls * D * F * M * (Nu+Nd)
    # Sparse object are currently stored as records, so here are their sizes:
    self.size_Vuu1 = dp + 4 * it
    self.size_Vdd1 = dp + 4 * it
    self.size_Vud1 = dp + 4 * it
    self.size_Ls1 = dp + 2 * it
    # memory required by the sparse objects, ALL determinants
    self.mem_Vss = (self.count_Vuu + self.count_Vdd) * self.size_Vuu1
    self.mem_Vxs = (self.count_Vud) * self.size_Vud1
    self.mem_Ls = (self.count_Ls) * self.size_Ls1

    if Print:
      self.printout_mem()

  def printout_mem(self, out=None):
    """
    Prints out a report for memory estimate.
    """

  # Tentative way to compute naive multithreading task sharing
  def task_div(self, D, th):
    """`Task divide-and-share':
    The division of D iterations into th threads ---
    to approximately account for imperfect task balance in OpenMP way.
    """
    inv_th = 1.0 / th
    return numpy.ceil(D * inv_th)

  def compute_step_cost(self, Print=False):
    # Placeholder: will be replaced by fancier stuff later
    # for more "symbolic" feel, or function that can give more actual
    # estimate of the operation cost.
    # For now these are merely
    LA = self.linalg
    mxm, mxv, vdot, mmtrace, tmmtrace = LA.mxm, LA.mxv, LA.vdot, LA.mmtrace, LA.tmmtrace
    (M, Nptot, Nu, Nd, F, D) = self.params_wlkr(0)
    try:
      th = self.num_threads
    except:
      th = 1

    d_fac = self.task_div(D, th)

    #self.cost_pre_Q = d_fac * F * mxm(M,N,N)
    #self.cost_Theta = d_fac * mxm(M,N,N)  # -- not considered for now

    self.ops_gf1_ovlp = d_fac * (mxm(Nu,Nu,M) + mxm(Nd,Nd,M)) # matmul of Psi^hc * Phi
    self.ops_gf1_ovlpinv = d_fac * (mxm(Nu,Nu,Nu) + mxm(Nd,Nd,Nd)) # the inverse of ovlp matrix
    self.ops_FB = d_fac * self.dens_Ls * F * (mmtrace(M,Nu) + mmtrace(M,Nd)) # the trace part
    self.ops_Elocal = d_fac * self.dens_Vss * (2*tmmtrace(M,M,Nu,Nu) + 2*tmmtrace(M,M,Nd,Nd) + tmmtrace(M,M,Nu,Nd)) # the trace part

    self.cost_gf1_ovlp = self.tpref_gf1_ovlp * self.ops_gf1_ovlp
    self.cost_gf1_ovlpinv = self.tpref_gf1_ovlpinv * self.ops_gf1_ovlpinv
    self.cost_FB = self.tpref_FB * self.ops_FB
    self.cost_Elocal = self.tpref_FB * self.ops_Elocal

    if Print:
      self.printout_compute()
