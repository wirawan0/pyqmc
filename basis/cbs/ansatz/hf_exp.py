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
      The exponential ansatz for SCF is attributed to Feller.

* CPL.286.243(1998)
  "Basis-set convergence in correlated calculations on Ne, N2, and H2O"
  Asger Halkier, Trygve Helgaker, Poul Jorgensen, Wim Klopper,
  Henrik Koch, Jeppe Olsen, Angela K. Wilson

  Oft-quoted in CCL mailing list for basis set extrapolation.
  But this paper discusses only the correlation energy extrapolation.

* JCP.108.154(1998)
  "An examination of intrinsic errors in electronic structure methods using
  the Environmental Molecular Sciences Laboratory computational results
  database and the Gaussian-2 set"
  Eq. 1


* CPL.302.437(1999)
  http://dx.doi.org/10.1016/S0009-2614(99)00179-7
  "Basis-set convergence of the energy in molecular Hartree-Fock calculations"
  Asger Halkier, Trygve Helgaker, Poul Jorgensen, Wim Klopper, Jeppe Olsen

  The original paper that proposed 2-point exponential-form
  CBS extrapolation of HF energies using c=1.63.

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
  def Guess_xy_old(self, x, y):
    imin = numpy.argmin(y)
    return (y[imin], 0, 0)
  def Guess_xy(self, x, y):
    """A better guess for exponential ansatz for CBS extrapolation.
    Naive guess with b = c = 0 can sometimes lead to pathological fit,
    so a better guess is needed.
    See routine Test_initial_guess for the rationale of the approximation.
    """
    from numpy import asarray, argsort, log, exp
    #xarr = asarray(x)
    #yarr = asarray(y)
    x = x[0]
    assert len(x) >= 3
    assert len(y) >= 3
    iorder = argsort(x)
    # choose the extremes in such a way to minimize the guess error
    (x1, x2, x3) = tuple(x[iorder[i]] for i in [0,-2,-1])
    (y1, y2, y3) = tuple(y[iorder[i]] for i in [0,-2,-1])
    D23 = y2 - y3
    D13 = y1 - y3
    E0 = y3
    c = -log(D23 / D13) / (x2 - x1)
    b = exp(c * x2) * (y2 - E0)
    return (E0, b, c)
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
  assert len(X) == len(E)
  rslt = DefObject().extrap(x=X, y=E)
  return rslt


""" --- Special two-value extrapolation for SCF energies ---

Example cases for Cr2 system (cc-pwCV[TQ5]Z-DK basis, RHF calculations):
see the analysis script, routine "examine_cbs_extrap_hf_20140728dset".

      E(X) = E_CBS + b * exp(-c*X)

                                                  -----------TQ5 fitting------------  -----------TQ fitting------------   -----------Q5 fitting------------
  rbond       E(TZ)       E(QZ)       E(5Z)       E_CBS       b           c           E_CBS       b           c           E_CBS       b           c
[ 1.5         0.69512178  0.69308933  0.69279279  0.69274213  0.76617248  1.92481764  0.69259408  0.33606784  1.63        0.69272053  0.25025876  1.63      ]
[ 1.55        0.70816632  0.7062469   0.70596556  0.70591724  0.71422788  1.92022707  0.70577919  0.31737792  1.63        0.70589701  0.23742841  1.63      ]
[ 1.6         0.72742565  0.72562349  0.72536297  0.72531894  0.69730714  1.93403249  0.72518436  0.29798789  1.63        0.72529948  0.21986646  1.63      ]
[ 1.6788      0.76579848  0.76423025  0.76397932  0.76393152  0.45569801  1.83250528  0.76384812  0.25930761  1.63        0.76391817  0.2117719   1.63      ]
[ 1.72        0.78832241  0.78686859  0.78661728  0.78656476  0.34028485  1.75526937  0.78651433  0.24038967  1.63        0.78655604  0.21208601  1.63      ]
[ 1.8         0.83472349  0.83346358  0.83321471  0.83315345  0.20370973  1.62186525  0.83315658  0.20832682  1.63        0.83315407  0.21002842  1.63      ]
[ 1.9         0.89488983  0.89380091  0.89355081  0.89347624  0.11666741  1.47106338  0.89353557  0.18005363  1.63        0.89348987  0.21107033  1.63      ]
[ 2.          0.95510086  0.95409784  0.95385247  0.95377301  0.09070589  1.40802156  0.95385343  0.16584983  1.63        0.95379268  0.20707089  1.63      ]
[ 2.4         1.1743032   1.1735533   1.17337398  1.17331762  0.0720804   1.43076928  1.17337057  0.12399692  1.63        1.17333028  0.15133372  1.63      ]
[ 3.          1.41304031  1.41267073  1.41256334  1.41251936  0.02123703  1.23594681  1.41258067  0.06111124  1.63        1.41253718  0.09062671  1.63      ]
[ 3.4         1.52227844  1.52177003  1.52176614  1.52176611  1.1496332e3 4.87458243  1.52164614  8.406640e-2 1.63        1.52176520  3.277315e-3 1.63      ]

The error E(CBS,TQ5) - E(CBS,TQ) seems to be small enough for our purpose.

"""

def extrap_cbs_hf_2pt(E, X=None, c=1.63):
  """Two-point CBS extrapolation for HF energy,
  according to

      E(X) = E_CBS + b * exp(-c*X)

  where only 2 values of X and E are available each.
  In psi4, for the Helgaker's two-parameter exponential ansatz,
  c is chosen to be 1.63.
  Let the X values be X1 and X2, where X2 = X1+1.

      E1 = E_CBS + b * exp(-c*X1)
      E2 = E_CBS + b * exp(-c*X1) * exp(-c)

      E1 - E2 = b * exp(-c*X1) (1 - exp(-c))

      b = (E1 - E2) / (exp(-c*X1) (1 - exp(-c)))

  Source: http://dx.doi.org/10.1016/S0009-2614(99)00179-7
  """
  from numpy import exp
  if X == None:
    X = numpy.arange(3, 3 + len(E), dtype=float)
  elif numpy.isscalar(X) and numpy.isreal(X):
    # ...if it is a numerical value
    X = numpy.arange(X, X + len(E), dtype=float)
  else:
    # ...otherwise, assume it is an array-like object
    X = numpy.array(X, dtype=float)

  assert len(X) == 2
  assert len(E) == len(X)

  X1, X2 = X
  assert abs(X1 - X2) == 1
  assert c > 0
  if X1 < X2:
    E1, E2 = E
  else:
    X2, X1 = X
    E2, E1 = E
  b = (E1 - E2) / (exp(-c*X1) * (1 - exp(-c)))
  E_CBS = E1 - b * exp(-c*X1)
  return E_CBS, b, c



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


def Test_initial_guess():
  """[20140329]

  The CBS extrap formula is:

      E(x) = E0 + b * exp(-c*x)

  Define:

      D12 = E(x1) - E(x2)
      D13 = E(x1) - E(x3)
      D23 = E(x2) - E(x3)

  Then:

      D12 = B exp(-c*x1) * ( 1 - exp(-c*(x2-x1)) )
      D13 = B exp(-c*x1) * ( 1 - exp(-c*(x3-x1)) )
      D23 = B exp(-c*x2) * ( 1 - exp(-c*(x3-x2)) )

      D23   1 - exp(-c*(x2-x3))
      --- = -------------------
      D13   1 - exp(-c*(x1-x3))

  which leads to an equation

      exp(-c*(x2-x1)) * ( 1 - (D12/D13) exp(-c*(x3-x2)) ) = 1 - (D12/D13)

  Now we do some approximation, noting the following:

  - if x1 < x2 < x3, then the dying exponential trend should lead to

      |D23| < |D12| < |D13|

  - furthermore c should be a positive constant

  The crudest approximation is this:

      D23   exp(-c*x2) - exp(-c*x3)
      --- = -----------------------  ~ exp(-c*(x2-x1))
      D13   exp(-c*x1) - exp(-c*x3)  ~

  This leads to a guess values of:
      E0 = min(E(:))
      c = -log( D23/D13 ) / (x2 - x1)
      b = exp(c*x2) * (E(x2) - E0)

  """
  from numpy import asarray
  from wpylib.text_tools import make_matrix as mtx
  Test_set_1_text = """
  # From Cr2 project, CBS extrap, r=1.6000
  # In the original algorithm < 20140328, this leads to unsuccessful fit.
  # x    Ebind (RHF)
    3    0.72742565
    4    0.72562349
    5    0.72536297
  """
  Test_set_1_params = [ # from gnuplot
    0.725318945019849, 0.697365649857968, 1.93406159078189
  ]
  dset_text = Test_set_1_text
  dset = asarray(mtx(dset_text))
  dset_x = dset[:,0]
  dset_y = dset[:,1]

  def make_guess(x, y):
    from numpy import asarray, argsort, log, exp
    #xarr = asarray(x)
    #yarr = asarray(y)
    assert len(x) >= 3
    assert len(y) >= 3
    iorder = argsort(x)
    # choose the extremes in such a way to minimize the guess error
    (x1, x2, x3) = tuple(x[iorder[i]] for i in [0,-2,-1])
    (y1, y2, y3) = tuple(y[iorder[i]] for i in [0,-2,-1])
    D23 = y2 - y3
    D13 = y1 - y3
    E0 = y3
    c = -log(D23 / D13) / (x2 - x1)
    b = exp(c * x2) * (y2 - E0)
    print E0, b, c
    return (E0, b, c)

  make_guess(dset_x, dset_y)

