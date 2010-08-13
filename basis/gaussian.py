# $Id: gaussian.py,v 1.3 2010-08-13 01:49:43 wirawan Exp $
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

----------
REFERENCES
----------

[SPHCART]
  Schlegel & Frisch, Int. J. Quant. Chem. 54, 83-87 (1995)


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
  # Map function type for packages other than GAMESS
  map_nwchem_func_type = { 'L': 'SP' }
  nfuncs_spher_map = {
    'S': 1, 'P': 3, 'D': 5, 'F': 7, 'G': 9, 'H': 11, 'I': 13,
    'L': 4,  # special (SP)
  }
  nfuncs_cart_map = {
    'S': 1, 'P': 3, 'D': 6, 'F': 10, 'G': 15, 'H': 21, 'I': 28,
    'L': 4,  # special (SP)
  }
  def __init__(self, species, srcfile=None, funcs=[], ispher=True):
    self.funcs = []
    self.species = species
    if funcs:
      self.funcs = funcs
      self.ispher = ispher
    elif srcfile:
      self.load(srcfile, ispher)
    # idea:
    # self.precision_control = GlobalPrecisionControl
    # and by changing local precision_control, we can alter the output
    # below.

  def load(self, srcfile, ispher=True, basisname=None):
    """Loads the basis function from a given source file.
    The data in the source file is essentially in GAMESS (US) format with
    minor modification to indicate the species and the end of the basis data
    (a line with only '.' in it)."""
    F = text_input(srcfile)
    try:
      F.seek_text(r'^(?i)spec: *%s$' % (self.species,))
    except StopIteration:
      raise ValueError, \
            "Basis functions not found for species `%s' in library file `%s'" \
            % (self.species, srcfile)
    if basisname == None:
      b = os.path.basename(srcfile)
      if b.endswith('.basis'): b = b[:-6]
      basisname = b
    self.basisname = b
    L = F.next_rec()
    funcs = []
    while L[0] != '.':
      # Reads a contracted basis
      # WARNING FIXME: L is still buggy. Still does not read the p coeffs correctly.
      func_typ = L[0].upper()
      func_gcount = int(L[1])
      func_desc = []
      for i in xrange(0, func_gcount):
        L = F.next_rec()
        func_desc.append((float(L[1]), float(L[2]))) # (exponent, coefficient) pair
      funcs.append((func_typ, func_desc))
      L = F.next_rec()
    self.funcs = funcs
    self.ispher = ispher

  def output_nwchem(self, specname=None, indent=0):
    """Outputs the basis definition in nwchem input format.
    The argument indent can be an integer (how many whitespaces added for
    indentation or a string."""
    specname = specname or self.species
    if isinstance(indent,int) and indent >= 0: indent = " " * indent

    return \
      "\n".join([
        "\n".join(["%s%s  %s" % (indent, specname, self.map_nwchem_func_type.get(typ,typ))] + [
            "%s  %17.8f %19.9f" % (indent, Exp, Coeff) \
              for (Exp,Coeff) in desc
        ]) \
          for (typ, desc) in self.funcs
      ])

  def output_psi3(self, specname=None, indent=0, basisname=None):
    """Outputs the basis definition in psi3 input format.
    The argument indent can be an integer (how many whitespaces added for
    indentation or a string."""
    specname = specname or self.species
    basisname = basisname or self.basisname
    if isinstance(indent,int) and indent >= 0: indent = " " * indent

    def list_sum(L):
      r = []
      for i in L:
        r += i
      return r

    return \
      "\n".join([
        '%s%s: "%s" = (' % (indent, specname, basisname),
        ] + list_sum([[
          '%s  (%s  ' % (indent, typ) + \
          ("\n%s      " % indent).join([
             "(%17.8f %19.9f)" % (Exp, Coeff) for (Exp,Coeff) in desc
          ]) + \
          ')'
        ] for (typ, desc) in self.funcs ]) + \
        [ '%s)' % (indent,), ]
      )

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

  def func_summary(self):
    """Gives a summary of the basis function, such as `9s8p5d1f'."""
    funcnames = 'SPDFGHI'
    count = dict([ (f,0) for f in funcnames ])
    for (typ, desc) in self.funcs:
      if typ == 'L':
        count['S'] += 1
        count['P'] += 1
      else:
        count[typ] += 1
    return "".join([ "%d%s" % (count[typ], typ.lower()) for typ in funcnames if count[typ] > 0 ])

  def func_count(self, ispher=None):
    """Gives the number of functions in this basis set.
    Each magnetic quantum number (for p, d, f, ... -type funcs) counts as
    a separate function in this formula."""
    if ispher == None: ispher = self.ispher
    if ispher:
      nfunc_map = self.nfuncs_spher_map
    else:
      nfunc_map = self.nfuncs_cart_map
    return sum([ nfunc_map[typ] for (typ, desc) in self.funcs ])



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

  def get(self, basis, species, fmt=None, specname=None, indent=0, ispher=True):
    if basis not in self:
      self[basis] = GTOBasisAtomList(basis)

    List = self[basis]
    if species not in List:
      List[species] = GTOBasis(species, srcfile=self.get_srcfile(basis), ispher=ispher)

    if fmt == None:
      return List[species] # raw, low-level data
    elif fmt == "nwchem":
      return List[species].output_nwchem(specname=specname, indent=indent)
    elif fmt == "gamess":
      return List[species].output_gamess(specname=specname, indent=indent)
    else:
      raise ValueError, "Invalid output format: " + fmt


# Global library: use this unless you need special customization, please.
LIB = GTOBasisLib()


def sph_norm(func):
  """This implements Eq.(8) in SPHCART paper, i.e. the normalization
  constant of a spherical Gaussian function.
  This will *NOT* do normalization for the 'L' type function"""
  (typ, desc) = func


