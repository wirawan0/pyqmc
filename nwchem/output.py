# $Id: output.py,v 1.2 2010-12-01 17:09:39 wirawan Exp $
#
# pyqmc.nwchem.output module
#
# Wirawan Purwanto
# Created: 20101028
#

"""
This module contains an output parser for nwchem code.

"""

import os
import os.path
import time

from wpylib.regexps import regex
from wpylib.iofmt.text_input import text_input
from wpylib.db.result_base import result_base
from pyqmc import PyqmcDataError

class nwchem_output(result_base):
  """Parser and structure for nwchem output text file.

  Fields from results_base:
    * filename_
    * absfilename_ -- full file name including absolute directory path.
    Here:
  Fields defined here:
    * info_code_version  nwchem version used in the calculation
    * info_mtime
    * title
  """

  class rx_:
    nwchem_version = regex(r'^\s*Northwest Computational Chemistry Package \(NWChem\)\s+([-0-9a-zA-Z_.]+)')
    # Obsoleted (not used)
    e_scf = regex(r'^\s*Total (SCF|DFT) energy *= *([-+eE0-9.]+)')
    e_ccsd = regex(r'^\s*CCSD total energy */ *hartree *= *([-+eE0-9.]+)')
    e_ccsd_t = regex(r'^\s*CCSD\(T\) total energy */ *hartree *= *([-+eE0-9.]+)')
    e_nucl = regex(r'^\s*Nuclear repulsion energy *= *([-+eE0-9.]+)')
    # Section headers
    hdr_input = regex(r'\s*NWChem Input Module\s*$')
    hdr_scf = regex(r'\s*NWChem SCF Module\s*$')
    hdr_mcscf = regex(r'\s*NWChem Direct MCSCF Module\s*$')
    underline = regex(r'\s*-+\s*$')
    # The following are not too foolproof, but nobody is so foolish to
    # use these as calculation title, isn't it?
    input_notitle = regex(r'^\s*Scaling coordinates for geometry ".*" by +[0-9.]+')
    scf_notitle = regex(r'^\s*ao basis += +".*"')

  def parse_file_(self, filename):
    """Extracts information from an nwchem output file (from its stdout).
    Right now, this parser is only good for single-point calculations
    (i.e. no multijob or geometry optimization at this point)."""

    Rx = self.rx_
    txtfile = text_input(filename,
                         skip_blank_lines=False,
                         comment_char='\0',
                         superize=True)

    self.clear()
    rslt = self
    self.txt_ = txtfile

    # This will also serve as an initial screening of the file
    try:
      L = txtfile.seek_text(Rx.nwchem_version.rx)
    except:
      raise PyqmcDataError, \
        "Cannot determine nwchem version in `%s'; perhaps it is not an nwchem output file" % filename
    rslt['info_code_version'] = (Rx.nwchem_version % L).group(1)
    rslt['info_mtime'] = time.localtime(os.path.getmtime(filename))

    search_patterns = [
      (regex(r'^\s*Total SCF energy *= *([-+eE0-9.]+)'),                   'E_SCF', float),
      (regex(r'^\s*Total MCSCF energy *= *([-+eE0-9.]+)'),                 'E_MCSCF', float),
      (regex(r'^\s*Total DFT energy *= *([-+eE0-9.]+)'),                   'E_DFT', float),
      (regex(r'^\s*CCSD total energy */ *hartree *= *([-+eE0-9.]+)'),      'E_CCSD', float),
      (regex(r'^\s*CCSD\(T\) total energy */ *hartree *= *([-+eE0-9.]+)'), 'E_CCSD_T', float),
      (regex(r'^\s*Nuclear repulsion energy *= *([-+eE0-9.]+)'),           'E_nuclear', float),
      (regex(r'^\s*Total MP2 energy\s+([-+eE0-9.]+)'),                     'E_MP2', float), # old MP2 module

      # The following are not that great.
      # We may want to use stricter scanner which takes into account the context
      # (which section) we are in
      #(regex(r'^\s*functions\s*=\s*([0-9]+)\s*$'),                         'nbasis', int),
    ]

    for L in txtfile:
      flds = L.split()
      if len(flds) == 0:
        continue
      elif self.detect_section_(Rx.hdr_input, L=L):
        # Input module output
        #print "L = ", L
        self.process_section_input_()
      elif self.detect_section_(Rx.hdr_scf, L=L):
        # SCF module output
        self.process_section_scf_()
      elif self.detect_section_(Rx.hdr_mcscf, L=L):
        # MCSCF module output
        self.process_section_mcscf_()
      else:
        for (pat, act, arg1) in search_patterns:
          if pat % L:
            if isinstance(act, basestring):
              rslt[act] = arg1(pat[1])
              break

    if 'nalpha_elec' in self:
      self['nelec_up'] = self.nalpha_elec
      self['nelec_dn'] = self.nbeta_elec

    return rslt


  def detect_section_(self, section_rx, F=None, L=None):
    """**Internal routine**

    Attempts to detect a section header text.
    Arguments:
    * section_rx = regex object for the section sought
    * F (optional) = the text_input object
    * L (optional) = the last text line read
    """
    if F == None: F = self.txt_
    if L == None:
      L1 = F.next()
    else:
      L1 = L

    if section_rx % L1:
      try:
        L2 = F.next()
      except:
        pass
      else:
        if self.rx_.underline % L2:
          return True
        else:
          F.file.push(L2)

    if L == None:
      F.file.push(L1)
    return False


  def skip_blank_lines_(self, num=None): #, min=None, max=None):
    txtfile = self.txt_
    if num != None:
      for i in xrange(num):
        L = txtfile.next()
        if L.strip() != "":
          raise PyqmcParseError, "Error skipping exactly %d blank lines." % num
    else:
      L = txtfile.next()
      while L.strip() == "": L = txtfile.next()
      txtfile.file.push(L)
    """ -- junked for now
    if min == None and max == None:
      c = 1
      while L.strip() == "":
        L = txtfile.next()
        c += 1
      txtfile.file.push(L)
      return c-1
    elif min == None:
      c = 0
      while c < max:
        if L.strip() != "":
          txtfile.file.push(L)
          break
        L = txtfile.next()
        c += 1
    """


  def process_section_input_(self):
    """**Internal routine**

    Process the output of 'input' module.
    """
    Rx = self.rx_
    txtfile = self.txt_

    # Look for job title, if any
    self.skip_blank_lines_(num=2)
    L = txtfile.next()
    #if not (Rx.input_notitle % L):
    if L.strip() != "" and not (Rx.input_notitle % L):
      self['title'] = L.strip()
      L = txtfile.next()  # strip underlining "-----"


  def process_section_scf_(self):
    """**Internal routine**

    Process the output of 'scf' module.
    """
    Rx = self.rx_
    txtfile = self.txt_

    search_patterns = [
      (regex(r'^\s*ao basis *= *"([^"]+)"'),                               'ao_basis', str),
      (regex(r'^\s*functions\s*=\s*([0-9]+)\s*$'),                         'nbasis', int),
      (regex(r'^\s*atoms\s*=\s*([0-9]+)\s*$'),                             'natoms', int),
      (regex(r'^\s*closed shells\s*=\s*([0-9]+)\s*$'),                     'nclosed_shells', int),
      (regex(r'^\s*open shells\s*=\s*([0-9]+)\s*$'),                       'nopen_shells', int),
      (regex(r'^\s*alpha electrons\s*=\s*([0-9]+)\s*$'),                   'nalpha_elec', int),
      (regex(r'^\s*beta electrons\s*=\s*([0-9]+)\s*$'),                    'nbeta_elec', int),
      (regex(r'^\s*charge\s*=\s*([0-9]+)\s*$'),                            'charge', int),
      (regex(r'^\s*wavefunction\s*=\s*([0-9]+)\s*$'),                      'scf_type', int),
    ]

    # Look for job title, if any
    self.skip_blank_lines_()
    L = txtfile.next()
    if not (Rx.scf_notitle % L):
      self['title_scf'] = L.strip()
      self.skip_blank_lines_()
    else:
      txtfile.file.push(L)

    for L in txtfile:
      flds = L.split()
      if len(flds) == 0:
        break
      else:
        for (pat, act, arg1) in search_patterns:
          if pat % L:
            if isinstance(act, basestring):
              self[act] = arg1(pat[1])
              break

    if 'nalpha_elec' not in self:
      # ROHF/RHF
      self['nalpha_elec'] = self.nclosed_shells + self.nopen_shells
      self['nbeta_elec'] = self.nclosed_shells
    else:
      # UHF
      pass


  def process_section_mcscf_(self):
    """**Internal routine**

    Process the output of 'mcscf' module.
    """
    Rx = self.rx_
    txtfile = self.txt_

    search_patterns = [
      (regex(r'^\s*Basis functions *: *([0-9]+)'),                         'nbasis_mcscf', int),
      (regex(r'^\s*Inactive shells\s*:\s*([0-9]+)\s*$'),                   'nfc', int),
      (regex(r'^\s*Active shells\s*:\s*([0-9]+)\s*$'),                     'nact_orb', int),
      (regex(r'^\s*Active electrons\s*:\s*([0-9]+)\s*$'),                  'nact_elec', int),
      (regex(r'^\s*Symmetry\s*:\s*([A-Za-z0-9]+)\s*$'),                    'symmetry_mcscf', str),
      (regex(r'^\s*Multiplicity\s*:\s*([0-9]+)\s*$'),                      'mult_mcscf', int),
    ]

    # Look for job title, if any
    self.skip_blank_lines_()
    L = txtfile.next()
    if not (Rx.underline % L):
      self['title_scf'] = L.strip()
      self.skip_blank_lines_()
    else:
      txtfile.file.push(L)

    for L in txtfile:
      flds = L.split()
      if len(flds) == 0:
        break
      else:
        for (pat, act, arg1) in search_patterns:
          if pat % L:
            if isinstance(act, basestring):
              self[act] = arg1(pat[1])
              break

    if 'nalpha_elec' not in self:
      self['nelec'] = 2 * self.nfc + self.nact_elec
      two_nbeta_elec = (self.nelec - self.mult_mcscf+1)
      if two_nbeta_elec % 2 != 0:
        raise PyqmcDataError, \
          "Invalid combination of multiplicity and num of active electrons?"
      self['nbeta_elec'] = two_nbeta_elec // 2
      self['nalpha_elec'] = self.nbeta_elec + self.mult_mcscf - 1


