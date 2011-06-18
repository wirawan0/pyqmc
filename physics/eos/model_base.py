# $Id: model_base.py,v 1.1 2011-06-18 02:52:11 wirawan Exp $
#
# pyqmc.physics.eos.model_base module
# Created: 20110414
# Wirawan Purwanto
#
# This module is part of PyQMC project.
#

"""
Module pyqmc.physics.eos.model_base

The base of all model equations.

SOME CONSIDERATIONS FOR EQUATIONS:

1. Always write the equations to allow scalar and array arguments.

2. Use numpy functions whenever possible.
"""

class model_base(object):
  """Base class of all EOS models."""
  def __init__(self, **_args):
    for (k,v) in _args.iteritems():
      setattr(self, k, v)

