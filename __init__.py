# $Id: __init__.py,v 1.3 2011-09-07 14:20:18 wirawan Exp $
#
# pyqmc module
# Created: ~200809
# Wirawan Purwanto
#
# Main pyqmc package. This will import everything.

#from pyqmc.multidet import *
#from pyqmc.matrices import *
#from pyqmc.walkers import *
pass

class PyqmcError(Exception):
  """PyQMC errors (unspecified in nature).
  This is a base error for more-specific PyQMC-related errors."""

class PyqmcDataError(PyqmcError):
  """Data related error (whether in-memory or on-disk)."""

class PyqmcParseError(PyqmcError):
  """Parser related error (e.g. for text files)."""

class PyqmcWarning(Warning):
  """PyQMC-related warnings."""

