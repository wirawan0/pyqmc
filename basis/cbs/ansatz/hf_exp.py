#
# pyqmc.basis.cbs.ansatz.hf_exp module
#
# Wirawan Purwanto
# Created: 20121011
#
"""
pyqmc.basis.cbs.ansatz.hf_exp module

Exponential form of ansatz used in CBS extrapolation,
mainly for Hartree-Fock energies.
See, e.g. JCP.106.9639(1997)
"Basis-set convergence of correlated calculations on water"
around Eq. (6), where the now-commonly-used way is described.

"""

import numpy
import wpylib
import wpylib.math.fitting

from wpylib.sugar import is_iterable

LITERATURE_TRACEBACK = """
Literature traceback:

* JCP.98.7059(1993)
  "The use of systematic sequences of wave functions for estimating the
  complete basis set, full configuration interaction limit in water"
  Under Sec. IV.

* JCP.99.1914(1993)
  "Benchmark calculations with correlated molecular wave functions.
  I. Multireference configuration interaction calculations for the second
  row diatomic hydrides"
  Sec. III. C., p. 1922
  Eq. 3

* JCP.106.9639(1997)
  "Basis-set convergence of correlated calculations on water"
  Eq. 6

      NOTE: This is the best example of paper using exp ansatz for HF
      CBS and 1/x^3 for correlation energy CBS.

* JCP.108.154(1998)
  "An examination of intrinsic errors in electronic structure methods using
  the Environmental Molecular Sciences Laboratory computational results
  database and the Gaussian-2 set"
  Eq. 1

"""

class hf_exp_fitting(object):
  """Exponential formula fitting for CBS extrapolation along
  correlation-consistent basis series,

     E(X) = E_CBS + b * exp(-c*X)

  This is mainly for 3-point extrapolation of Hartree-Fock energies.
  """

  NOTES_1 = """
  This was imported from my Ca+4H2_Hstorage_iter.py script
  CVS rev 1.108.
  Original class name: CBS_exp_fitting.

  Fitting coefficients:
  * C[0] = Total energy at CBS limit
  * C[1] = linear coefficient of exp
  * C[2] = in exponent

  As with SHO and AHO3, leastsq does not perform well.
  The best way is to fit with `fmin' simplex method.
  """
  dim = 1  # a function with 1-D domain
  fit_opts = dict(xtol=1e-8, maxfun=100000, maxiter=10000, disp=0)
  debug = 0
  def __call__(self, C, x):
    return C[0] + C[1] * numpy.exp(-numpy.abs(C[2])*x[0])
  def Guess_xy(self, x, y):
    imin = numpy.argmin(y)
    return (y[imin], 0, 0)
  def extrap(self, x, y):
    """Main function for basis extrapolation.

    Input:
    * `x' values (as an array) are the cc-pVxZ basis cardinal number
    * `y' values are the corresponding total (HF) energy.
    """
    if len(x.shape) == 1:
      # fix common mistake: make it 2-D
      x = x.reshape((1, x.shape[0]))
    self.last_fit = wpylib.math.fitting.fit_func(
                      Funct=self,
                      x=x, y=y, method='fmin',
                      debug=self.debug,
                      outfmt=0, # yield full result
                      opts=self.fit_opts,
                    )
    return self.last_fit['xopt']
  def extrap345(self, e3, e4, e5):
    """Simple function to do CBS extrapolation in TZ-QZ-5Z
    sequence."""
    x = numpy.array([[3., 4., 5.]])
    y = numpy.array([e3, e4, e5])
    return self.extrap(x, y)

#CBS_exp = CBS_exp_fitting()

# Convenience functions:

def DefObject():
  """Returns a default, globally-defined instance
  of the extrapolation/fitting class defined in this module."""
  global Def_hf_exp_fitting
  try:
    return Def_hf_exp_fitting
  except:
    pass
  Def_hf_exp_fitting = hf_exp_fitting()
  return Def_hf_exp_fitting


def extrap_cbs_hf(E, X=None):
  """Performs CBS many-point exponential extrapolation in X.
  This is originally intended for HF energy only.

  * If X is None, then the series is assumed to start from
    3 (TZ, QZ, ...).
  * If X is a scalar number, then it is taken to be the starting
    point of the series.
  * Otherwise it must be an array-like object containing the
    array values corresponding to the energies in E.
  """
  if X == None:
    X = numpy.arange(3, 3 + len(E), dtype=float)
  elif numpy.isscalar(X) and numpy.isreal(X):
    # ...if it is a numerical value
    X = numpy.arange(X, X + len(E), dtype=float)
  else:
    # ...otherwise, assume it is an array-like object
    X = numpy.array(X, dtype=float)
  rslt = DefObject().extrap(x=X, y=E)
  return rslt



# Tests below (for my sanity)

def Test_1():
  """
  Source:
  Ca+4H2, cc-pCV[TQ5]Z basis

   X  E(UHF)              File
   3  -681.070870436283   $GAFQMC/nwchem/Ca+4H2/Z2.3000/cc-pCVTZ/UHF/dHH0.7682/Ca+4H2_UHF_Z2.3000.dHH0.7682.cc-pCVTZ.06.out
   4  -681.073942347902   $GAFQMC/nwchem/Ca+4H2/Z2.3000/cc-pCVQZ/UMP2/dHH0.7682/Ca+4H2_UMP2_Z2.3000.dHH0.7682.cc-pCVQZ.06.out
   5  -681.074733235173   $GAFQMC/nwchem/Ca+4H2/Z2.3000/cc-pCV5Z/UHF/dHH0.7682/Ca+4H2_UHF_Z2.3000.dHH0.7682.cc-pCV5Z.06.out


  # Note:  Ca+4H2_UHF_Z2.3000.dHH0.7682.cc-pCVQZ.06.out contains slightly higher UHF
  # energy (-681.072975269138) because of near-linear-dependency.

  >>> hf_exp.Test_1()
  array([ -6.81075007e+02,   2.42420813e-01,   1.35689988e+00])

  # More detailed values:
  # E_cbs = -681.07500745524953
  # b     = 0.24242081344723787
  # c     = 1.3568998808614841
  """
  E = [-681.070870436283, -681.073942347902, -681.074733235173]
  x = 3
  return extrap_cbs_hf(E,x)


def Test_2():
  """
  Source:
  Ca+4H2, cc-pCV[TQ5]Z basis
  Similar to above, but result from numbers in my GROFF notebook.

   X  E(UHF)
   3  -681.070870436
   4  -681.073942348
   5  -681.074733235


  Text from Ca+4H2-hydrogen-storage.groff: (l.5769 of CVS rev 1.90):
           Z    d(H-H)     E(TZ)           E(QZ)           E(5Z)           E(CBS)
  Ca4H2    2.3  0.7682     -681.070870436  -681.073942348  -681.074733235  -681.075007454812  b=0.24242114379679 c=1.35690034749527

  >>> hf_exp.Test_2()
  array([ -6.81075007e+02,   2.42421146e-01,   1.35690035e+00])

  # More detailed values:
  # E_cbs = -681.07500745480684
  # b     = 0.24242114585492372
  # c     = 1.3569003507370354

  # Note: this verifies OK against the GROFF notebook results above.
  """
  E = [-681.070870436,  -681.073942348,  -681.074733235]
  x = 3
  return extrap_cbs_hf(E,x)
