#
# pyqmc.basis.cbs.ansatz.corr_x3 module
#
# Wirawan Purwanto
# Created: 20121011
#
"""
pyqmc.basis.cbs.ansatz.corr_x3 module

1/x**3 form of ansatz used in CBS extrapolation,
mainly for correlation energies.
See, e.g. JCP.106.9639(1997)
"Basis-set convergence of correlated calculations on water"
around Eq. (7), where the now-commonly-used way is described.

"""

import numpy
import wpylib
import wpylib.math.fitting

from wpylib.sugar import is_iterable
from wpylib.math.stats.errorbar import errorbar
from wpylib.math.fitting.linear import linregr2d_SZ
from wpylib.math.fitting import fit_result


LITERATURE_TRACEBACK = """
Literature traceback:

* JCP.106.9639(1997)
  "Basis-set convergence of correlated calculations on water"
  Eqs. 6 & 7

      NOTE: This is the best example of paper using exp ansatz for HF
      CBS and 1/x^3 for correlation energy CBS.

* T. Muller, "Basis sets, accuracy, and calibration in quantum chemistry",
  in "Computational Nanoscience: Do It Yourself!", NIC Series,
  Vol. 31, edited by J. Grotendorst, S. Blugel, and D.Marx
  (John von Neumann Institute for Computing, Julich, 2006) pp. 1943.

"""

is_numeric_scalar = lambda X: numpy.isscalar(X) and numpy.isreal(X)

class corr_x3_fitting(object):
  """Exponential formula fitting for correlation energy CBS
  extrapolation along the correlation-consistent basis series,

     Ecorr(X) = Ecorr(CBS) + B/(X**3)

  where X is the cardinal number in cc-pVxZ series.
  """

  NOTES_1 = """
  This was imported from my Ca+4H2_Hstorage_iter.py script
  CVS rev 1.108.
  Original class name: none (methods was from class cbx_class).

  Fitting coefficients:
  * C[0] = Total energy at CBS limit
  * C[1] = linear coefficient of 1/X**3
  """
  dim = 1  # a function with 1-D domain
  fit_opts = dict(xtol=1e-8, maxfun=100000, maxiter=10000, disp=0)
  debug = 0
  def __call__(self, C, x):
    return C[0] + C[1] * (x[0]**(-3))
  def Guess_xy(self, x, y):
    imin = numpy.argmin(y)
    return (y[imin], 0)
  def extrap(self, x, y):
    """Main function for basis extrapolation.

    Input:
    * `x' values (as an array) are the cc-pVxZ basis cardinal number
    * `y' values are the corresponding correlation energy.
      Note that y can contain errorbar (i.e. values given as errorbar
      objects).
    """
    x = numpy.array(x, copy=False).flatten()
    y = numpy.array(y, copy=False).flatten()
    invx3 = x**(-3)

    if isinstance(y[0], errorbar):
      # special handling
      if len(y) == 2:
        rslt = extrap_2pt(y[0], x[0], y[1], x[1])
        # energy and slope already in rslt
        self.last_fit = fit_result(
          xopt=rslt,
        )
      else:
        # more general 3+ point linear extrapolation
        yy = [ y1.mean() for y1 in y ]
        dy = [ y1.error() for y1 in y ]
        #print invx3
        #print yy
        #print dy
        r = linregr2d_SZ(invx3, yy, dy)
        #print r
        rslt = (errorbar(r['a'], r['sigma_a']),
                errorbar(r['b'], r['sigma_b']),
               )
        #print rslt
        self.last_fit = r
        r['xopt'] = rslt
    else:
      # the old way for scalars
      rslt = numpy.polyfit(x=invx3, y=y, deg=1, full=True)
      keys = ('xopt', 'residuals', 'rank', 'singular_values', 'rcond')
      self.last_fit = fit_result(dict(zip(keys,rslt)))

    return self.last_fit['xopt']

# Convenience functions:

def DefObject():
  """Returns a default, globally-defined instance
  of the extrapolation/fitting class defined in this module."""
  global Def_corr_x3_fitting
  try:
    return Def_corr_x3_fitting
  except:
    pass
  Def_corr_x3_fitting = corr_x3_fitting()
  return Def_corr_x3_fitting


def extrap_cbs_corr(E, X=None):
  """Performs CBS many-point extrapolation of correlation
  energy in X.
  This is originally intended for correlation energy only.

  * If X is None, then the series is assumed to start from
    3 (TZ, QZ, ...).
  * If X is a scalar number, then it is taken to be the starting
    point of the series.
  * Otherwise it must be an array-like object containing the
    array values corresponding to the energies in E.
  """
  if X == None:
    X = numpy.arange(3, 3 + len(E), dtype=float)
  elif is_numeric_scalar(X):
    # ...if it is a numerical value
    X = numpy.arange(X, X + len(E), dtype=float)
  else:
    # ...otherwise, assume it is an array-like object
    X = numpy.array(X, dtype=float)
  rslt = DefObject().extrap(x=X, y=E)
  return rslt


def extrap_2pt(E1, x1, E2, x2):
  """Performs CBS two-point linear extrapolation in 1/x**3.
  This tool should properly account the error bars."""
  dx = 1.0 / x1**3 - 1.0 / x2**3
  dE = E1 - E2
# # The simple formula is:
# slope_CBS = dE / dx
# E_CBS = E2 - slope_CBS / x2**3
  #     = E2 - (E1-E2) / ( 1.0 / x1**3 - 1.0 / x2**3 ) / x2**3
  #    := E2 - (E1-E2) * fac
  #     = E2 * (1 + fac) - E1 * fac
  # The alternate formula is: -- this should factor in the errorbar
  fac = 1.0 / dx / x2**3
  E_CBS = E2 * (1 + fac) - E1 * fac
  return (E_CBS, slope_CBS)



# Even more-convenient functions--these are cross-cutting concerns:

def extrap_cbs(E, E_HF, X=None, X_HF=None, debug=0):
  """
  Does the 'standard' CBS extrapolation using separate methods
  for HF and correlation energies:

  * Feller's exp() ansatz for the HF energy CBS extrapolation
  * 1/X**3 ansatz for correlation energy CBS extrapolation

  This requires the energies being calculated on a
  correlation-consistent basis family.

  As an example of an actual application, please consult
  J. Chem. Phys. 106, 9639 (1997)
  (http://link.aip.org/link/doi/10.1063/1.473863).


  NOTE: E and E_HF can be a dict, in which case X
  and X_HF arguments should be ignored.
  The dict keys are the corresponding `X' argument.
  """
  from pyqmc.basis.cbs.ansatz.hf_exp import extrap_cbs_hf

  if debug >= 20:
    print "E      = ", E
    print "E_HF   = ", E_HF
    print "X      = ", X
    print "X_HF   = ", X_HF

  if isinstance(E, dict):
    e = E
  else:
    if X == None:
      X = numpy.arange(3, 3 + len(E), dtype=float)
    elif is_numeric_scalar(X):
      X = numpy.arange(X, X + len(E), dtype=float)
    # recast E to a dict
    e = dict(zip(X, E))

  if isinstance(E_HF, dict):
    e_hf = E_HF
  else:
    if X_HF == None:
      X_HF = numpy.arange(3, 3 + len(E_HF), dtype=float)
    elif is_numeric_scalar(X_HF):
      X = numpy.arange(X_HF, X_HF + len(E_HF), dtype=float)
    # recast E to a dict
    e_hf = dict(zip(X_HF, E_HF))

  e_corr = dict([ (x1, e[x1] - e_hf[x1]) for x1 in e.keys() ])
  if debug >= 10:
    print "e_corr =", e_corr

  x_hf = sorted(e_hf.keys())
  x_corr = sorted(e_corr.keys())
  if len(x_hf) < 3:
    if debug >= 0:
      print "extrap_cbs: Warning: HF has fewer than three data points; taking", \
        "the largest X=%s as the CBS limit" % (x_hf[-1])
    cbs_hf = [ e_hf[x_hf[-1]] ]
  else:
    cbs_hf = extrap_cbs_hf([ e_hf[x1] for x1 in x_hf ], x_hf)

  cbs_corr = extrap_cbs_corr([ e_corr[x1] for x1 in x_corr ], x_corr)

  if debug >= 10:
    print "cbs_hf =", cbs_hf
    print "cbs_corr =", cbs_corr
  return (cbs_hf[0] + cbs_corr[0], cbs_corr[1], cbs_hf, cbs_corr)



# Tests below (for my sanity)

def Test_1():
  """
  Source:
  Ca+4H2, cc-pCV[TQ5]Z basis, Z=2.3 cc-pCVxZ basis

   X  E(UHF)             E(QMC/UHF)      File
   3  -681.070870436283  -681.61620(61)  $GAFQMC/nwchem/Ca+4H2/Z2.3000/cc-pCVTZ/UHF/dHH0.7682/Ca+4H2_UHF_Z2.3000.dHH0.7682.cc-pCVTZ.06.out
   4  -681.073942347902  -681.83561(74)  $GAFQMC/nwchem/Ca+4H2/Z2.3000/cc-pCVQZ/UMP2/dHH0.7682/Ca+4H2_UMP2_Z2.3000.dHH0.7682.cc-pCVQZ.06.out
   5  -681.074733235173  -681.9189(12)   $GAFQMC/nwchem/Ca+4H2/Z2.3000/cc-pCV5Z/UHF/dHH0.7682/Ca+4H2_UHF_Z2.3000.dHH0.7682.cc-pCV5Z.06.out

  QMC numbers from my GROFF notebook

  # Note:  Ca+4H2_UHF_Z2.3000.dHH0.7682.cc-pCVQZ.06.out contains slightly higher UHF
  # energy (-681.072975269138) because of near-linear-dependency.

  >>> corr_x3.Test_1()
  Correlated energies:
    -0.54533(61)
    -0.76167(74)
    -0.8442(12)
  Total energy @ CBS =  -681.9976(10)
  (errorbar(-0.922559144089,0.00103522799244,'-0.9226(10)'),
   errorbar(10.1930295272,0.0365555298231,'10.193(37)'))

  """
  from wpylib.math.stats.errorbar import errorbar

  eb = errorbar.create_str
  E_HF = [-681.070870436283, -681.073942347902, -681.074733235173]
  x = 3
  E_cbs_HF = -681.07500745524953
  #cbs_HF = extrap_cbs_hf(E,x) - no need; get value from
  # pyqmc.basis.cbs.ansatz.hf_exp.Test_1() output!
  E_QMC = [ eb("-681.61620(61)"), eb("-681.83561(74)"), eb("-681.9189(12)") ]
  E_corr = [ (e2-e1) for (e2,e1) in zip(E_QMC, E_HF) ]
  print "Correlated energies: "
  for e1 in E_corr:
    print " ", e1
  cbs_corr = extrap_cbs_corr(E_corr, [3,4,5])
  print "Total energy @ CBS = ", cbs_corr[0] + E_cbs_HF
  return cbs_corr


def Test_2():
  """
  Similar to Test_1 but exercises the cross-cutting extrap_cbs().

  >>> corr_x3.Test_2()
  Correlated energies:
    -0.54533(61)
    -0.76167(74)
    -0.8442(12)
  E      =  [errorbar(-681.6162,0.00061,'-681.61620(61)'), errorbar(-681.83561,0.00074,'-681.83561(74)'), errorbar(-681.9189,0.0012,'-681.9189(12)')]
  E_HF   =  [-681.07087043628303, -681.07394234790195, -681.07473323517297]
  X      =  [3, 4, 5]
  X_HF   =  [3, 4, 5]
  e_corr = {3: errorbar(-0.545329563717,0.00061,'-0.54533(61)'), 4: errorbar(-0.761667652098,0.00074,'-0.76167(74)'), 5: errorbar(-0.844166764827,0.0012,'-0.8442(12)')}
  cbs_hf = [ -6.81075007e+02   2.42420813e-01   1.35689988e+00]
  cbs_corr = (errorbar(-0.922559144089,0.00103522799244,'-0.9226(10)'), errorbar(10.1930295272,0.0365555298231,'10.193(37)'))
  Total energy @ CBS =  -681.9976(10)
  (errorbar(-681.997566599,0.00103522799244,'-681.9976(10)'),
   errorbar(10.1930295272,0.0365555298231,'10.193(37)'),
   array([ -6.81075007e+02,   2.42420813e-01,   1.35689988e+00]),
   (errorbar(-0.922559144089,0.00103522799244,'-0.9226(10)'),
    errorbar(10.1930295272,0.0365555298231,'10.193(37)')))

  """
  from wpylib.math.stats.errorbar import errorbar

  eb = errorbar.create_str
  E_HF = [-681.070870436283, -681.073942347902, -681.074733235173]
  E_QMC = [ eb("-681.61620(61)"), eb("-681.83561(74)"), eb("-681.9189(12)") ]
  E_corr = [ (e2-e1) for (e2,e1) in zip(E_QMC, E_HF) ]
  print "Correlated energies: "
  for e1 in E_corr:
    print " ", e1
  cbs = extrap_cbs(E_QMC, E_HF, [3,4,5], [3,4,5], debug=20)
  print "Total energy @ CBS = ", cbs[0]
  return cbs

