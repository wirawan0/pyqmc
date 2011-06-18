# $Id: pec_morse.py,v 1.1 2011-06-18 02:52:11 wirawan Exp $
#
# pyqmc.physics.pec_morse module
# Created: 20110414
# Wirawan Purwanto
#
# This module is part of PyQMC project.
#

"""
Module pyqmc.physics.eos.pec_morse

This module contains details of 1-D, Morse-like
potential energy curve models.
"""

import numpy
from pyqmc.physics.eos.model_base import model_base
from numpy import exp

class morse(model_base):
  """Standard Morse potential energy equation:

      morse(x) = E0 + De * (1 - exp(-a * (x - re)))**2

  Required parameters:

  E0, De, a, re

  References:
  1.  http://en.wikipedia.org/wiki/Morse_potential
  2.  P. M. Morse, "Diatomic molecules according to the wave mechanics.
      II. Vibrational levels."
      Phys. Rev. 1929, 34, 57-64.
      doi:10.1103/PhysRev.34.57
  """
  _params_ = ("E0", "De", "a", "re")
  def __call__(self, x):
    return self.E0 + self.De * (1 - exp(-self.a * (x - self.re)))**2


class morse2(model_base):
  """Morse potential energy equation with a twist:

      morse2(x) = E0 + 0.5 * k / a**2 * (1 - exp(-a * (x - re)))**2

  Required parameters:

  E0, k, a, re

  This twist is devised to obtain the spring constant (k) more easily.
  It would also allow estimating uncertainties of k better.
  This identifies De with:

     De = 0.5 * k / a**2
  """
  _params_ = ("E0", "De", "a", "re")
  def __call__(self, x):
    return self.E0 + 0.5 * self.k * ( (1 - exp(-self.a * (x - self.re))) / self.a )**2

