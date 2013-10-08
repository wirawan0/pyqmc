# -*- python -*-
#
# pyqmc.utils.linalg_cost module
#
# Wirawan Purwanto
# Created: 20131004
#

"""
pyqmc.utils.linalg_cost

Cost estimator and analyzer for linear algebra operations.

"""

class linalg_cost_base(object):
  """A base class for all other linear-algebra cost estimators.
  """

class linalg_cost_ops(linalg_cost_base):
  """Simple estimator based only on the number of FMA operations.
  """
  def mxm(self, M,N,K):
    """Matrix - matrix multiply:
        A(M,K) . B(K,N)
    """
    return M*K*N
  def mxv(self, M,N):
    """Matrix - vector multiply."""
    return M*N
  def vdot(self, M):
    """Vector - vector dot product."""
    return M
  def mmtrace(self, M, N):
    """Trace of a matrix-matrix product:
        Tr(A(M,N) . B(N,M))
    """
    return M*N
  def tmmtrace(self, M, N, P, Q):
    """Trace of a tensor-matrix-matrix product:
        Tr(V(M,N,P,Q) . A(M,P) . B(N,Q))
    """
    return M*N*P*Q



# TODO: estimate or provide a way to [empirically] account for
# mem cache effects
