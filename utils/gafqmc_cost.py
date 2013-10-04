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

from pyqmc.utils import cost


class gafqmc_cost_estimator(cost.qmc_cost_estimator):
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

  # Turns out, empirically Vss and Vxs densities are the same


  def compute_mem_cost_sparse1(self, Print=False):
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
    (M, Nptot, F, D) = (self.nbasis, self.nptot, self.nflds, self.npsitdet)
    try:
      Nu = self.nup
      Nd = self.ndn
    except:
      print >> sys.stderr, "WARNING: cannot find nup/ndn, making a ballpark estimate"
      Nu = (self.nptot + 1) // 2
      Nd = self.nptot - Nu
      print >> sys.stderr, "WARNING: estimated nup, ndn =", nup, ndn
    (dpc, dp, it) = (self.size_complex, self.size_real, self.size_int)

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

