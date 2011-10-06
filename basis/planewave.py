# $Id: planewave.py,v 1.1 2011-10-06 19:41:50 wirawan Exp $
#
# pyqmc.basis.planewave module
# Created: 20100929
# Wirawan Purwanto
#
# This is part of PyQMC project.
#
# Module for planewave basis
#

"""
Module for planewave basis.

"""

import numpy

class grid3d(object):
  """A three-dimensional grid for plane wave (or real-space) basis."""
  def __init__(self, LL):
    """
    LL is a three-dimensional vector containing the number of grid points
    in each dimension.
    """

