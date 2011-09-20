# $Id: gaussian.py,v 1.8 2011-09-20 22:01:46 wirawan Exp $
#
# pyqmc.basis.gaussian module
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
import wpylib.file.file_utils
from wpylib.file.file_utils import path_search
import wpylib.iofmt.text_input
from wpylib.iofmt.text_input import text_input
from wpylib.sugar import ifelse, list_join

# TODO: GlobalPrecisionControl object which controls the precision
#

class GTOBasis(object):
  """
  Gaussian-type orbital basis.
  Attributes:
  . fmt_exp
  . fmt_coeff
  . funcs: [ (angmom, [(exp1,coeff1),(exp2,coeff2),...]), ... ]
  . species
  . ispher
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
  fmt_exp = "%17.8f"    # format for exponent printout
  fmt_coeff = "%19.9f"  # format for coefficient printout
  def __init__(self, species, srcfile=None, funcs=[], ispher=True):
    self.funcs = []
    self.species = species
    self.ispher = ispher
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
    (a line with only '.' in it).

    An example of valid field is like this:

      Spec: H
      S   3
        1     33.8650000  0.0254938
        2      5.0947900  0.1903730
        3      1.1587900  0.8521610
      S   1
        1      0.3258400  1.0000000
      S   1
        1      0.1027410  1.0000000
      S   1
        1      0.0360000  1.0000000
      P   1
        1      0.7500000  1.0000000
      .

    """
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
      if func_typ == 'L' or func_typ == 'SP':
        # For now, the 'L' or 'SP' function is loaded into two terms,
        # one s and one p.
        # This corresponds to the way nwchem does this thing.
        func_desc_s = []
        func_desc_p = []
        for i in xrange(0, func_gcount):
          L = F.next_rec()
          func_desc_s.append((float(L[1]), float(L[2]))) # (exponent, coefficient) pair
          func_desc_p.append((float(L[1]), float(L[3]))) # (exponent, coefficient) pair
        funcs.append(('S', func_desc_s))
        funcs.append(('P', func_desc_p))
      else:
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

    FMT = "%%s  %s %s" % (self.fmt_exp, self.fmt_coeff)
    return \
      "\n".join([
        "\n".join(["%s%s  %s" % (indent, specname, self.map_nwchem_func_type.get(typ,typ))] + [
            FMT % (indent, Exp, Coeff) \
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

    FMT = "(%s %s)" % (self.fmt_exp, self.fmt_coeff)
    return \
      "\n".join([
        '%s%s: "%s" = (' % (indent, specname, basisname),
        ] + list_join(*[[
          '%s  (%s  ' % (indent, typ) + \
          ("\n%s      " % indent).join([
             FMT % (Exp, Coeff) for (Exp,Coeff) in desc
          ]) + \
          ')'
        ] for (typ, desc) in self.funcs ]) + \
        [ '%s)' % (indent,), ]
      )

  def output_mpqc(self, specname=None, indent=0, basisname=None):
    """Outputs the basis definition in mpqc input format.
    The argument indent can be an integer (how many whitespaces added for
    indentation or a string."""
    specname = specname or self.species
    basisname = basisname or self.basisname
    if isinstance(indent,int) and indent >= 0: indent = " " * indent

    FMT = "%s %s" % (self.fmt_exp, self.fmt_coeff)
    return \
      "\n".join([
        '%s%s: "%s": [' % (indent, specname, basisname),
        ] + list_join(*[[
          '%s  (type: [(am=%s%s)]\n' % (indent, typ.lower(),
                                       ifelse(self.ispher and self.nfuncs_spher_map[typ] >= 5, " puream=1", "")) + \
          '%s   {exp coef:0} = {' % (indent,) + \
          ("\n%s      " % indent).join([""] + [
             FMT % (Exp, Coeff) for (Exp,Coeff) in desc
          ]) + \
          '\n%s  })' % (indent,)
        ] for (typ, desc) in self.funcs ]) + \
        [ '%s]' % (indent,), ]
      )

  def output_gamess(self, specname=None, indent=0):
    """Outputs the basis definition in gamess(US) input format.
    The argument indent can be an integer (how many whitespaces added for
    indentation or a string.
    Note that due to GAMESS way of specifying basis functions, the atom
    name is not included."""
    specname = specname or self.species
    if isinstance(indent,int) and indent >= 0: indent = " " * indent

    FMT = "%%s %%3d %s %s" % (self.fmt_exp, self.fmt_coeff)
    return \
      "\n".join([
        "\n".join(["%s%s  %d" % (indent, typ, len(desc))] + [
            FMT % (indent, idx+1, Exp, Coeff) \
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
  basis function name.

  Attributes:
  . name: basis name
  . ispher (if available, determines default "ispher" for newly created
    GTOBasis objects)
  . fmt_exp,
  . fmt_coeff: if defined, used to override class GTOBasis's defaults when
    creating a new object.
  """
  def __init__(self, name, ispher=True):
    self.name = name
    self.ispher = ispher
  def get(self, species, srcfile=None, ispher=None):
    """Gets a GTOBasis class for a given species.
    If not existing, will create a new one."""
    if ispher == None: ispher = self.ispher
    if species not in self:
      self[species] = GTOBasis(species, srcfile=srcfile, ispher=ispher)
      if hasattr(self, "fmt_exp"): self[species].fmt_exp = self.fmt_exp
      if hasattr(self, "fmt_coeff"): self[species].fmt_coeff = self.fmt_coeff
    return self[species]



class GTOBasisLib(dict):
  """Library of basis functions, categorized by basis name (e.g. cc-pVTZ) and
  species (e.g. Ca).
  The basis functions are loaded on as-needed basis and cached in memory."""
  def __init__(self):
    import pyqmc
    self.path = [ os.path.join(os.path.dirname(pyqmc.__file__), "libs/basis") ]
    self.fname_map = {
      '6-31++G**': '6-31++Gxx',
      '6-311++G**': '6-311++Gxx',
    }
    pass

  def get_srcfile(self, basis):
    return path_search(self.path, self.fname_map.get(basis, basis) + ".basis")

  def get(self, basis, species, fmt=None, specname=None, indent=0, ispher=True):
    """FIXME: Definition of ispher here is only used when creating a new,
    non-existent object."""
    if basis not in self:
      self[basis] = GTOBasisAtomList(basis)

    List = self[basis]
    if species not in List:
      # will create a new basis object if needed:
      srcfile = self.get_srcfile(basis)
      if srcfile == None:
        raise ValueError, "Cannot find library file for basis: %s" % basis
      List.get(species, srcfile=srcfile, ispher=ispher)
      #List[species] = GTOBasis(species, srcfile=self.get_srcfile(basis), ispher=ispher)

    if fmt == None:
      return List[species] # raw, low-level data
    elif fmt == "nwchem":
      return List[species].output_nwchem(specname=specname, indent=indent)
    elif fmt == "gamess":
      return List[species].output_gamess(specname=specname, indent=indent)
    elif fmt == "psi3":
      return List[species].output_psi3(specname=specname, indent=indent)
    elif fmt == "mpqc":
      return List[species].output_mpqc(specname=specname, indent=indent)
    else:
      raise ValueError, "Invalid output format: " + fmt


# Global library: use this unless you need special customization, please.
LIB = GTOBasisLib()


def sph_norm(func):
  """This implements Eq.(8) in SPHCART paper, i.e. the normalization
  constant of a spherical Gaussian function.
  This will *NOT* do normalization for the 'L' type function"""
  (typ, desc) = func


def convert_BSE_basis(fname, outfname):
  """Converts basis information obtained from BSE to internal format.
  This subroutine expect GAMESS-US format."""

  from time import strftime
  from pyqmc.physics.atoms import get as get_atom

  In = text_input(fname)
  In.skip_blank_lines = False

  if outfname == str:
    class outwriter:
      def __init__(self):
        self.r = []
      def write(self,s):
        self.r.append(s)
      def writelines(self,s):
        self.r += s
      def __str__(self):
        return "".join(self.r)
    out = outwriter()
  else:
    out = open(outfname, "w")

  out.writelines(["# Converted automatically using pyqmc.basis.gaussian.convert_BSE_basis\n",
                  "# Source file: %s\n" % fname,
                  "# Date: %s\n" % strftime("%Y-%m-%d %H:%M"),
                 ])
  for I in In:
    if I.startswith("!"):
      I = "#" + I[1:]
    elif I.startswith("$DATA"):
      break # go to the next stage
    out.writelines([I, "\n"])

  ended = 0
  for spec in In:
    atom = get_atom(spec)
    out.writelines(["Spec: ", atom.symb, "\n"])
    for bas in In:
      if bas.strip() == "":
        out.write(".\n\n")
        break
      elif bas.startswith("$END"):
        out.write(".\n\n")
        ended = 1
        break
      else:
        out.writelines([bas, "\n"])
    if ended:
      break

  if (outfname == str):
    return str(out)
