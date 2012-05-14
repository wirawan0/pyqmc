#
# pyqmc.abinit.output module
#
# Wirawan Purwanto
# Created: 20120217
#

"""
This module contains an output parser for abinit output.

"""

import numpy
import os
import os.path
import sys
import time
import weakref

from wpylib.params.params_struct import Struct as struct
from wpylib.regexps import regex
from wpylib.iofmt.text_input import text_input
from wpylib.iofmt.text_output import text_output
from wpylib.db.result_base import result_base
from pyqmc import PyqmcDataError
from warnings import warn

# see http://www.abinit.org/documentation/helpfiles/for-v6.10/input_variables/vardev.html#istwfk
istwfk_special_kpt = [
  (2, (0.0, 0.0, 0.0)),
  (3, (0.5, 0.0, 0.0)),
  (4, (0.0, 0.0, 0.5)),
  (5, (0.5, 0.0, 0.5)),
  (6, (0.0, 0.5, 0.0)),
  (7, (0.5, 0.5, 0.0)),
  (8, (0.0, 0.5, 0.5)),
  (9, (0.5, 0.5, 0.5)),
]

DEBUG_LEVEL = 0

def is_special_kpt(kpt, istwfk):
  for (istwfk1, kpt1) in istwfk_special_kpt:
    if (istwfk == 0 or istwfk == istwfk1) and (kpt == kpt1):
      return True
  return False

def get_num_pw(npw, kpt, istwfk):
  if is_special_kpt(kpt, istwfk):
    if kpt == (0,0,0):
      return npw * 2 + 1
    else:
      return npw * 2
  else:
    return npw


class abinit_dataset(result_base):
  """A single dataset of ABINIT calculation.

  Important results:
  - scf: True or False.
  - scf_nsteps
  - scf_converged
  - Etotal
  - nkpt
  - wfk[spin][ikpt] = eigensolution of a particular spin sector (`up' or `dn')
    and k-point (0, 1, ...nkpt-1)

  IMPORTANT: For consistency with python-style processing, all indices are
  converted to python-based (e.g. the k-point numbering starts from zero
  instead of one)."""
  class rx_:
    # This sneaky info can be used to determine # of planewaves:
    kpt_npw = regex(r'P newkpt: treating\s+(?P<nband>[0-9]+) bands with npw=\s*(?P<npw>[0-9]+) for ikpt=\s*(?P<ikpt>[0-9]+)')
    scf_begin = regex(r'^\s*iter\s+Etot\(hartree\)\s+deltaE\(h\)') # marker of SCF block
    scf_line1 = regex(r'^\s*ETOT\s*[0-9\*]+\s+(?P<Etot>[-+eE0-9.]+)\s+(?P<Ediff>[-+eE0-9.]+)')
    scf_convg1 = regex(r'^\s*At SCF step\s+(?P<numscf>[0-9]+)\s*,\s*etot is converged')
    E_begin = regex(r'^\s*Components of total free energy')
    E_end = regex(r'^\s*-+\s*$')
    E_parts_list = [
      (regex(r'^\s*Kinetic energy\s*=\s*([-+eE0-9.]+)'), 'E_kinetic', float),
      (regex(r'^\s*Hartree energy\s*=\s*([-+eE0-9.]+)'), 'E_hartree', float),
      (regex(r'^\s*Ewald energy\s*=\s*([-+eE0-9.]+)'), 'E_ewald', float),
      (regex(r'^\s*PspCore energy\s*=\s*([-+eE0-9.]+)'), 'E_pspcore', float),
      (regex(r'^\s*Loc\.?\s+psp\.?\s+energy\s*=\s*([-+eE0-9.]+)'), 'E_psploc', float),
      (regex(r'^\s*NL\s+psp\.?\s+energy\s*=\s*([-+eE0-9.]+)'), 'E_pspnonloc', float),
      (regex(r'^\s*>+ Internal E\s*=\s*([-+eE0-9.]+)'), 'E_internal', float),
      (regex(r'^\s*-kT\*entropy\s*=\s*([-+eE0-9.]+)'), 'E_smear', float),
      (regex(r'^\s*>+ Etotal\s*=\s*([-+eE0-9.]+)'), 'Etotal', float),
    ]
    eigen_begin = regex(r'^\s*Eigenvalues \(hartree\) for nkpt=\s*(?P<nkpt>[0-9]+)\s*k points(?:, SPIN (?P<spin>[A-Za-z]+))?:')
    eigen_kpt1 = regex(r'^\s*kpt#\s*(?P<ikpt>[0-9]+),\s*nband=\s*(?P<nband>[0-9]+), wtk=\s*(?P<wtk>[.0-9]+), kpt=\s*(?P<kx>[-+.0-9]+)\s+(?P<ky>[-+.0-9]+)\s+(?P<kz>[-+.0-9]+)')
    eigen_kpt_stop = regex(r'^\s*prteigrs\s*: .*do not print more k-points')
    densph_begin = regex(r'^\s*Atom\s+Sphere radius\s+Integrated_up_density\s+Integrated_dn_density\s+Total[^\s]+\s+Diff')
    dataset_end = regex(r'^\s*==\s*(?:DATASET\s+[0-9]+|END DATASET\(S\))\s+=+')

  class dt_:
    scf_cycle = numpy.dtype([('Etot', float), ('deltaE', float)])

  def parse_dataset_results_(self, baserec, F, index):
    """Parses a dataset: the `result' part.
    The """

    scf_cycle = []
    rx = self.rx_
    dt = self.dt_
    try:
      dbg = DEBUG_FD
    except:
      dbg = text_output(None) # sys.stdout, flush=True)
    self['MyErrors'] = 0
    self['dataset_index'] = index
    self['parent_'] = weakref.ref(baserec)

    # We can get the following from the outvars section of the output preamble
    def getoutvar(kwd, default, mapfunc=None):
      try:
        r = getattr(baserec.outvars, kwd + str(index))
      except AttributeError:
        r = getattr(baserec.outvars, kwd, default)
      try:
        return map(mapfunc, r)
      except:
        if mapfunc == None:
          return r
        else:
          try:
            return mapfunc(r)
          except:
            return r
    #ngfft = getattr(baserec.outvars, "ngfft%d" % index, getattr(baserec.outvars, "ngfft", None))
    self['ngfft'] = tuple(getoutvar("ngfft", None, mapfunc=int))
    istwfk = getoutvar("istwfk", None, mapfunc=int)
    kpt_npw = {}

    for L in F:
      dbg("L:  %s\n" % L.rstrip())
      if rx.dataset_end % L:
        dbg("** end dataset detected **\n")
        F.file.push(L) # put back the text data to the file
        break
      elif rx.kpt_npw % L:
        # record the number of planewaves (w/o istwfk correction)
        kpt_npw[int(rx.kpt_npw['ikpt'])-1] = int(rx.kpt_npw['npw'])
      elif rx.scf_begin % L:
        # Extracts the SCF cycle data plus whether it converges
        dbg("** SCF section **\n")
        for L2 in F:
          #dbg("L2: %s\n" % L2.rstrip())
          if not (rx.scf_line1 % L2):
            break
          scf_cycle.append( (float(rx.scf_line1['Etot']), float(rx.scf_line1['Ediff'])) )
        self['scf_data'] = numpy.array(scf_cycle, dtype=dt.scf_cycle)
        self['scf'] = True

        for L2 in F:
          if L2.strip() != "": break

        if rx.scf_convg1 % L2:
          self['scf_nsteps'] = rx.scf_convg1['numscf']
          self['scf_converged'] = True

      elif rx.E_begin % L:
        dbg("** Energy section **\n")
        for L2 in F:
          L2 = L2.rstrip()
          dbg("L2: %s\n" % L2)
          if rx.E_end % L2:
            break
          elif len(L2) == 0:
            continue
          else:
            for (pat, act, arg1) in rx.E_parts_list:
              if pat % L2:
                if isinstance(act, basestring):
                  self[act] = arg1(pat[1])
                  break

      elif rx.eigen_begin % L:
        # Begins reading info of eigensolutions here: eigenvals, occ, ...
        self['nkpt'] = int(rx.eigen_begin['nkpt'])
        if (istwfk == None): istwfk = (0,) * self['nkpt']
        spin = rx.eigen_begin['spin']
        dbg("** Eigenvector sector: spin = %s **\n" % spin)
        if spin == None:
          self['udet'] = False
          spins = ('up',)
          wfk = result_base(up=[])
          wfk['dn'] = wfk['up']
        elif str(spin).lower() == 'up':
          self['udet'] = True
          spins = ('up', 'dn')
          wfk = result_base(up=[], dn=[])
        else:
          raise PyqmcDataError, "Error: Unknown spin type detected"
        self['wfk'] = wfk

        for (ispin, s) in enumerate(spins):
          wfk_s = wfk[s]
          for k in xrange(self['nkpt']):
            L2 = F.next()
            dbg("  spin %s kpt %d: %s\n" % (s, k, L2.strip()))
            if not (rx.eigen_kpt1 % L2):
              if (rx.eigen_kpt_stop % L2):
                warn("Cannot parse all kpt data; reason: `%s'" % L2)
                self['MyErrors'] += 1
                break
              else:
                raise PyqmcDataError, "Expected `kpt#' line, got `%s'" % L2
            e = rx.eigen_kpt1
            # CAVEAT: The results here have a very limited precision.
            wfk_s_k = result_base(
              ikpt = int(e['ikpt'])-1, # convert to 0-based
              nband = int(e['nband']),
              wtk = float(e['wtk']),
              kpt = (float(e['kx']), float(e['ky']), float(e['kz']),),
            )
            npwbasis = get_num_pw(kpt_npw[k], wfk_s_k.kpt, istwfk[k])
            if DEBUG_LEVEL >= 10:
              print (kpt_npw[k], wfk_s_k.kpt, istwfk[k]), "->", npwbasis
            wfk_s_k['npwbasis'] = npwbasis
            wfk_s_k['istwfk'] = istwfk[k]
            if k != wfk_s_k['ikpt']:
              warn("Unexpected kpt index in Abinit output: given %d, expecting %d" \
                   % (wfk_s_k['ikpt'], k))
              self['MyErrors'] += 1
              # In case of discrepancy, we proceed, but at your risk
            # Reads in eigenvalues
            egnvals = []
            mbands = wfk_s_k['nband']
            for L2 in F:
              dbg("  egn: " + L2)
              m = map(float, L2.split())
              egnvals += m
              mbands -= len(m)
              if mbands <= 0: break
            wfk_s_k['egnval'] = numpy.array(egnvals)

            # Reads in occupancy (if any)
            L2 = F.next()
            if not (L2.strip().startswith('occupation numbers')):
              # Maybe this one does not have occupancy; push back and continue
              # FIXME: Must get the occ from the end of the file.
              # For now we don't give the 'occ' field.
              F.file.push(L2)
            else:
              occ = []
              mbands = wfk_s_k['nband']
              for L2 in F:
                m = map(float, L2.split())
                occ += m
                mbands -= len(m)
                if mbands <= 0: break
              wfk_s_k['occ'] = numpy.array(occ)

            wfk_s.append(wfk_s_k)

          if ispin+1 < len(spins):
            assert rx.eigen_begin % F.next()

      elif rx.densph_begin % L:
        # spin density within a spherical boundary: for local spin moments
        self['local_spin_moment'] = \
          F.read_items((0, int, 'iatom'), (1, float, 'sphere_radius'),
                       (2, float, 'up_spin_dens'), (3, float, 'dn_spin_dens'),
                       (4, float, 'tot_spin_dens'), (5, float, 'diff_spin_dens'),
                       end_line_match=r'Note: Diff')

      else:
        pass

  @property
  def parent(self):
    p = self.parent_()
    if p == None:
      raise ReferenceError, "Parent output data structure was already deleted."
    else:
      return p


class abinit_output(result_base):
  """Parser and structure for abinit main output text file
  (the ab_out text file, not the stdout content).

  Fields from result_base:
    * filename_
    * absfilename_ -- full file name including absolute directory path.
    Here:
  Fields defined here:
    * info_code_version  abinit version used in the calculation
    * info_mtime
  """

  class rx_:
    abinit_version = regex(r'^[. ]Version\s+(?P<version>[^\s]+)\s+of ABINIT')
    dataset_begin = regex(r'^\s*==\s*DATASET\s+(?P<dataset>[0-9]+)\s+=+')
    outvars_begin = regex(r'-outvars: echo values of preprocessed input variables')

  def parse_file_(self, filename):
    """Extracts information from an abinit text output file.
    Right now, this parser is only good for single-point calculations
    (i.e. no multijob or geometry optimization at this point)."""

    rx = self.rx_

    self.clear()

    txtfile = text_input(filename, comment_char='\0', skip_blank_lines=False, superize=1)
    # This will also serve as an initial screening of the file
    try:
      L = txtfile.seek_text(rx.abinit_version.rx)
    except:
      raise RuntimeError, \
        "Cannot determine ABINIT version in `%s'; perhaps it is not an ABINIT output file" % filename
    self['info_code_version'] = (rx.abinit_version % L).group('version')
    self['info_mtime'] = time.localtime(os.path.getmtime(filename))

    dataset = {}
    self['dataset'] = dataset
    for L in txtfile:
      if rx.dataset_begin % L:
        dset = int(rx.dataset_begin['dataset'])
        dataset[dset] = abinit_dataset()
        dataset[dset].parse_dataset_results_(self, txtfile, index=dset)
      elif rx.outvars_begin % L:
        self.parse_outvars_(txtfile)

    return self

  def parse_outvars_(self, txtfile):
    """Parses the outvars section at the beginning of abinit output.
    """
    Vars = struct()
    self.outvars = Vars
    x = struct()
    x.vname = None
    x.vval = []
    def vflush(x):
      if x.vname != None:
        setattr(Vars, x.vname, tuple(x.vval))
      x.vname = None
      x.vval = []
    for L in txtfile:
      F = L.split()
      if len(F) == 0:
        break
      if F[1] == ':':
        # some comments show up like this:
        #   outvar1 : prtvol=0, do not print more k-points.
        #   prtocc : prtvol=0, do not print more k-points.
        # we ignore those
        continue
      if F[0] == 'P':
        del F[0]
      if F[0][0].isalpha():
        vflush(x)
        x.vname = F[0]
        x.vval = F[1:]
      else:
        x.vval += F
    vflush(x)


