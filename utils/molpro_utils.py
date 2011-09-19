# $Id: molpro_utils.py,v 1.1 2011-09-19 20:07:37 wirawan Exp $

"""
pyqmc.utils.molpro_utils
Collection of MOLPRO-related utilities.

Wirawan Purwanto
"""

import numpy
import os
import sys
import re

from pyqmc import PyqmcDataError
from pyqmc.physics import atoms
import pyqmc.basis.gaussian as gbs

from wpylib.iofmt.text_input import text_input
from wpylib.regexps import regex as rx
from wpylib.params import struct


class molpro_text_input(text_input):
  """Molpro-style text_input."""
  comma_Rx = re.compile(r'(?:\s+,\s*|\s*,\s+)')
  def __init__(self, *_argsl, **_argsd):
    super(molpro_text_input,self).__init__(*_argsl, **_argsd)
    self.comment_char = "!"
    self.set_next_proc(self.next_line_molpro)
  def next_line_molpro(self):
    """Reads a complete text line from molpro.
    If it continues to the next line, we slurp the next line as well."""
    rslt = []
    while True:
      LL = self.next_line()
      L2 = self.comma_Rx.sub(",", LL)
      rslt.append(L2)
      if not L2.endswith(","):
        break
    return "".join(rslt)


class molpro_basis_file(gbs.GTOBasisAtomList):
  """PyQMC representation of MOLPRO Basis file."""
  def read(self, infile, clobber=True, floatize=True):
    """Parses basis text in MOLPRO format.
    Options:
    - clobber: if True, will erase all basis definitions
      before reading in the data from the input file.
    - floatize: if True, converts the number to float.
      Otherwise will preserve them as strings.

    NOTES:
    In http://theochem.weizmann.ac.il/web/papers/group12.html ,
    the authors seem to use ``%.13f'' format to display their
    MOLPRO-formatted tables.
    """
    if clobber: self.clear()
    inp = molpro_text_input(infile)
    self._readstate = struct(angmom=None, exps=None, contracted=None,
                             infile=infile, inp=inp)
    state = self._readstate
    if floatize:
      state.Float = float
    else:
      state.Float = lambda x: x

    for L in inp:
      F = L.split(',')
      state.F = F
      what = F[0].lower()
      if what == "c":
        #print "1:", L
        if state.angmom == None:
          raise PyqmcDataError, \
            "%s:%d: Illegal contraction: there was no exponents defined yet." \
            % (infile, inp.lineno)
        try:
          (c_begin, c_end) = tuple([ int(d)-1 for d in F[1].split('.') ])
        except:
          raise PyqmcDataError, \
            "%s:%d: Unparseable contraction range string: %s" \
            % (infile, inp.lineno, F[1])
        self._add_contraction(c_begin, c_end, F[2:])
        state.contracted = True
      elif len(what) == 1 and what in "spdfghi":
        #print "2:", L
        self._flush_exponents()  # flush uncontracted exponents if any
        spec = F[1]
        state.angmom = what.upper()
        state.spec = atoms.get(spec).symb  # makes sure symbol is valid
        state.exps = F[2:] # [ float(f1) for f1 in F[2:] ]
        #print ">>", state.exps
        state.contracted = False
      else:
        raise PyqmcDataError, \
          "%s:%d: Invalid line detected: %s" \
          % (infile, inp.lineno, L)
    self._flush_exponents()  # flush remaining uncontracted exponents if any
    inp.close()
    del self._readstate

  def _add_contraction(self, c_begin, c_end, coeffs):
    """** Internal routine for molpro basis reader. **
    Add contraction based on a given set coefficients and contraction
    range limits."""
    state = self._readstate
    #print "add: ", c_begin, c_end, coeffs
    #print "add::", state.exps
    if c_begin < 0 or c_begin >= len(state.exps) \
       or c_end < c_begin or c_end >= len(state.exps):
      raise PyqmcDataError, \
        "%s:%d: Invalid contraction range limits: %s" \
        % (state.infile, state.inp.lineno, state.F[1])
    if len(coeffs) != c_end - c_begin + 1:
      raise PyqmcDataError, \
        "%s:%d: Coefficients array length (%d) disagrees with contraction limits (%s)" \
        % (state.infile, state.inp.lineno, len(coeffs), state.F[1])
    basisobj = self.get(state.spec)
    func_desc = [ (state.Float(e), state.Float(c))
                    for (e,c) in zip(state.exps[c_begin:c_end+1], coeffs) ]
    #print state.exps[c_begin-1:c_end-1]
    #print coeffs
    #print zip(state.exps[c_begin-1:c_end-1], coeffs)
    basisobj.funcs.append((state.angmom, func_desc))
    state.contracted = True

  def _flush_exponents(self):
    """** Internal routine for molpro basis reader. **
    Flush list of exponents as uncontracted functions.
    This is the default unless there are contractions defined."""
    state = self._readstate
    if state.angmom == None or state.contracted:
      return
    basisobj = self.get(state.spec)
    basisobj.funcs += [ (state.angmom, [ (state.Float(f),1.0) ] )
                          for f in state.exps ]
    state.angmom = None
    state.exps = None
    state.contracted = None
    



