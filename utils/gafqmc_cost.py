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


class gafqmc_cost_estimator(cost_base):
  """Cost estimator for GAFQMC calculation.

  Names of precomputed objects:
  * Vijkl = four-indexed two-body operator
  * Vss = product of two trial wfn orbitals (same spin) with Vijkl
  * Vxs = product of two trial wfn orbitals (opposite spins) with Vijkl

  """

  # (SPARSE) PRECOMPUTED MATRICES

  # Default sparse matrix density (ballpark estimate)
  # These values MAY NOT BE CORRECT for your particular calculation!
  # These are ok for GTO-basis calculations without frozen core:
  dens_Vss = 0.5
  dens_Vxs = 0.5
  dens_rho = 0.7

  def compute_mem_cost(self, Print=False):
    """
    Estimate a calculation's MEMORY cost based on the given input sizes.

    Required input:
    - nbasis
    - nflds
    - nwlkmax
    - nptot, nup, ndn
    - npsitdet

    Objects:

    Output:
    - mem_Vss =
    - mem_V
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
    (dpc, dp) = (self.size_complex, self.size_real)

    self.wlk_size = Nptot * M * dpc
    self.mem_wlk = self.wlk_size * nwlkmax_proc
    self.mem_Lvec = self.get_hsop_dim * F * dp   # so far it is double precision
    """ UNDER CONSTRUCTION
    self.mem_Vss =
    """

    if Print: raise NotImplementedError

  def printout_mem(self, out=None):
    """
    Prints out a report for memory estimate.
    """

