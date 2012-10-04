#!/usr/bin/python
# $Id: __init__.py,v 1.1 2009-06-12 15:06:50 wirawan Exp $
#
# pyqmc.utils module
#
# Wirawan Purwanto
# Created: 20090329
#

import os
import os.path


def _file_search(Dir, Names):
  """Look for one of the files listed in `Names` in subdirectory `Dir`.
  The first one found will be returned.
  This is useful in case the file of interested may be compressed or
  have different names (from different generation of tools/programs used)."""
  for N in Names:
    fn = os.path.join(Dir, N)
    if os.path.isfile(fn):
      return fn
  return None

