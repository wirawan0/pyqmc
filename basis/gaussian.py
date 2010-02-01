# $Id: gaussian.py,v 1.1 2010-02-01 21:29:16 wirawan Exp $
#
# pyqmc.basis module
# Created: 20100201
# Wirawan Purwanto
#
# This is part of PyQMC project.
#
# Module for GTO basis
#
# Todo:
# * Right now I output 1 more decimal digit than the BSE output.
#   But precision control for the output routines need to be
#   added.

"""
Module for GTO basis.
"""

import os.path
import re
import pyqmc.utils.file_utils
from pyqmc.utils.file_utils import path_search
import pyqmc.utils.text_input
from pyqmc.utils.text_input import text_input

# TODO: GlobalPrecisionControl object which controls the precision
#

class GTOBasis(object):
  """
  Gaussian-type orbital basis.
  """
  def __init__(self, species, srcfile=None, funcs=[]):
    self.funcs = []
    self.species = species
    if funcs:
      self.funcs = funcs
    elif srcfile:
      self.load(srcfile)
    # idea:
    # self.precision_control = GlobalPrecisionControl
    # and by changing local precision_control, we can alter the output
    # below.

  def load(self, srcfile):
    F = text_input(srcfile)
    try:
      F.seek_text(r'^(?i)spec: *%s$' % (self.species,))
    except StopIteration:
      raise ValueError, \
            "Basis functions not found for species `%s' in library file `%s'" \
            % (self.species, srcfile)
    L = F.next_rec()
    funcs = []
    while L[0] != '.':
      # Reads a contracted basis
      func_typ = L[0]
      func_gcount = int(L[1])
      func_desc = []
      for i in xrange(0, func_gcount):
        L = F.next_rec()
        func_desc.append((float(L[1]), float(L[2]))) # (exponent, coefficient) pair
      funcs.append((func_typ, func_desc))
      L = F.next_rec()
    self.funcs = funcs

  def output_nwchem(self, specname=None, indent=0):
    """Outputs the basis definition in nwchem input format.
    The argument indent can be an integer (how many whitespaces added for
    indentation or a string."""
    specname = specname or self.species
    if isinstance(indent,int) and indent >= 0: indent = " " * indent

    return \
      "\n".join([
        "\n".join(["%s%s  %s" % (indent, specname, typ)] + [
            "%s  %17.8f %19.9f" % (indent, Exp, Coeff) \
              for (Exp,Coeff) in desc
        ]) \
          for (typ, desc) in self.funcs
      ])


  def output_gamess(self, specname=None, indent=0):
    """Outputs the basis definition in gamess(US) input format.
    The argument indent can be an integer (how many whitespaces added for
    indentation or a string.
    Note that due to GAMESS way of specifying basis functions, the atom
    name is not included."""
    specname = specname or self.species
    if isinstance(indent,int) and indent >= 0: indent = " " * indent

    return \
      "\n".join([
        "\n".join(["%s%s  %d" % (indent, typ, len(desc))] + [
            "%s %3d %17.8f %19.9f" % (indent, idx+1, Exp, Coeff) \
              for (idx, (Exp,Coeff)) in zip(xrange(len(desc)), desc)
        ]) \
          for (typ, desc) in self.funcs
      ])


class GTOBasisAtomList(dict):
  """List of basis functions (GTOBasis) for many species, for the same
  basis function name."""
  def __init__(self, name):
    self.name = name


class GTOBasisLib(dict):
  """Library of basis functions, categorized by basis name (e.g. cc-pVTZ) and
  species (e.g. Ca).
  The basis functions are loaded on as-needed basis and cached in memory."""
  def __init__(self):
    import pyqmc
    self.path = [ os.path.dirname(pyqmc.__file__) + "/libs/basis" ]
    self.fname_map = {
      '6-311++G**': '6-311++Gxx',
    }
    pass

  def get_srcfile(self, basis):
    return path_search(self.path, self.fname_map.get(basis, basis) + ".basis")

  def get(self, basis, species, fmt=None, specname=None, indent=0):
    if basis not in self:
      self[basis] = GTOBasisAtomList(basis)

    List = self[basis]
    if species not in List:
      List[species] = GTOBasis(species, srcfile=self.get_srcfile(basis))

    if fmt == None:
      return List[species]
    elif fmt == "nwchem":
      return List[species].output_nwchem(specname=specname, indent=indent)
    elif fmt == "gamess":
      return List[species].output_gamess(specname=specname, indent=indent)
    else:
      raise ValueError, "Invalid output format: " + fmt


# Global library: use this unless you need special customization, please.
LIB = GTOBasisLib()
