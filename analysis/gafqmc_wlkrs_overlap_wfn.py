#!/usr/bin/python
#
# pyqmc.analysis.gafqmc_wlkrs_overlap_wfn module
#
# Wirawan Purwanto
# Created: 20150402
#


"""
pyqmc.analysis.gafqmc_wlkrs_overlap_wfn

[20150402]
An analysis module to examine overlap between GAFQMC walker(s) and
another wave function.


"""

import numpy
import os
import sys
from wpylib.iofmt.text_output import text_output
from pyqmc.matrices import slaterdet
from pyqmc.data.gafqmc_walkers import wlk_gafqmc
from pyqmc.results.gafqmc_info import gafqmc_info


class gafqmc_wlkrs_overlap_wfn(object):
  """
  Main analysis class for examining overlap between GAFQMC walker(s) and
  another wave function.
  """
  _chkfile_pattern_openmp = "%(wlkr_dir)s/W_old_run_000"
  _chkfile_pattern_legacy = "%(wlkr_dir)s/W_old_run_%(rank)05d"
  _chkfile_pattern_new = "%(wlkr_dir)s/gafqmc-%(rank)05d"
  _chkfile_pattern_blk_suffix = ".b%(measblk)010d"

  # defaults:
  chkfile_pattern = _chkfile_pattern_legacy
  chkfile_pattern_blk_suffix = _chkfile_pattern_blk_suffix

  def init(self, info, wlkr_dir, chkfile_pattern=None, logfile=sys.stdout):
    self._init_metadata(info)
    self.wlkr_dir = wlkr_dir
    if chkfile_pattern is not None:
      self.chkfile_pattern = chkfile_pattern
    self.logfile = text_output(logfile)

  def _init_metadata(self, info):
    if isinstance(info, basestring):
      self.infofile = info
      info = gafqmc_info(info)
    self.nup = info.nup
    self.ndn = info.ndn
    self.udet = info.udet
    self.nbasis = info.nbasis
    self.H0 = info.H0
    try:
      self.num_tasks = info.num_tasks
    except:
      self.num_tasks = 1
    try:
      self.info_file = info.info_file
    except:
      pass

  def load_all_walkers(self, measblk=None, verbose=10):
    """Loads a set of checkpoint files (either the only one set
    [with blk=None], or that belonging to a particular snapshot of
    the walker population.
    """
    O = self.logfile
    O("load_all_walkers: Reading %d walker files; wlkr_dir=%s, measblk=%s\n" \
      % (self.num_tasks, self.wlkr_dir, measblk))
    self.chkdata = Chk = wlk_gafqmc()
    Chk.init_metadata(nbasis=self.nbasis, nup=self.nup, ndn=self.ndn,
                      udet=self.udet)

    args = dict(wlkr_dir=self.wlkr_dir)
    if measblk is not None:
      chkfile_pattern = self.chkfile_pattern + self.chkfile_pattern_blk_suffix
      args['measblk'] = measblk
    else:
      chkfile_pattern = self.chkfile_pattern

    # data_all is the result data obtained from reading all the walkers,
    # including the QMC walkers
    Chk.data_all = Chk.read_all(self.num_tasks,
                                self.chkfile_pattern,
                                fname_args=args,
                                verbose=verbose)
    Chk.wlkrs_all = Chk.data_all.wlkrs  # quick access
    self.wlkrs_all = Chk.data_all.wlkrs  # quick access
    O("load_all_walkers: Total number of dets read in: %d\n" \
      % Chk.wlkrs_all.ndets)
    return Chk

  def dets_assign_amplitudes(self):
    dets_assign_amplitudes(self.wlkrs_all)

  def dets_orthonormalize(self):
    self.wlkrs_all.normalize_dets()



def dets_compute_impfn(wlkrs, PsiT):
  """Given a multideterminant wave function from QMC checkpoint
  (i.e. a QMC random walker population), recalculate
  the important function <PsiT|phi> for each determinant |phi>
  in the population.

  This step is needed when the importance function was not saved into the
  checkpoint file.
  """
  from pyqmc.matrices.slaterdet import mdet_ovlp, Det, MultiDet
  D1 = Det()
  D1.ampl = 1.0
  MD = MultiDet()
  MD.dets = [D1]
  for D in wlkrs.dets:
    D1.up, D1.dn = D.up, D.dn
    ovlp_PsiT_phi = mdet_ovlp(PsiT, MD)
    D.impfn = ovlp_PsiT_phi


def dets_assign_amplitudes(wlkrs):
  """Given a multideterminant wave function from QMC checkpoint,
  assign the amplitude.
  """
  if not hasattr(wlkrs.dets[0], 'impfn'):
    print("Warning: checkpoint data does not have impfn field!")
    for D in wlkrs.dets:
      D.ampl = D.wtwlkr
  else:
    for D in wlkrs.dets:
      D.ampl = D.wtwlkr / D.impfn


def orthonormalize_dets(wlkrs):
  """Orthonormalizes the orbitals in each Slater determinant in the
  multideterminant series `wlkrs`.

  This is an experimental routine to be used in conjunction to a
  gafqmc_wlkrs_overlap_wfn object.

  Note: for walker population in gafqmc_wlkrs_overlap_wfn object,
  this functionality has been superseded by the .dets_normalize() method of
  that class, provided that the amplitude (.ampl field of each walker)
  has been assigned by calling the .dets_assign_amplitudes() method.

  Thus this standalone function is somewhat obsolete; only call when
  situation really warrants it.
  """
  from wpylib.math.linalg import modgs
  from numpy import product
  for D in wlkrs.dets:
    """Attributes of the determinant:
    - up    (alpha)
    - dn    (beta)
    - ampl  (supposed to be: wtwlkr / impfn * det_norm )
    From checkpoint:
    - impfn (if available)
    - wtwlkr
    - det_norm, up_det_norm, dn_det_norm
    - up_orb_norms, dn_orb_norms

    """
    # treat everything like udet case only
    # should not be harmful other then space/time efficiency waste
    up_det, D.up_orb_norms = modgs(D.up)
    dn_det, D.dn_orb_norms = modgs(D.dn)

    D.up_det_norm = product(D.up_orb_norms)
    D.dn_det_norm = product(D.dn_orb_norms)
    D.det_norm = D.up_det_norm * D.dn_det_norm
    D.ampl = getattr(D, 'ampl', 1.0) * getattr(D, 'wtwlkr', 1.0) * D.det_norm / getattr(D, 'impfn', 1.0)
    D.up = up_det
    D.dn = dn_det


def load_fort70_wfn(chkdata, fort70file, verbose=10):
  """Loads a (multideterminant) wave function from fort.70 file.
  This is necessary to get the wave function in the orthogonal basis.
  Another alternative would be to reconstitute the X basis-orthogonalization
  matrix from fort.70 and transform the wave function orbitals manually to the
  orthogonal basis.
  """
  from numpy import zeros, eye, asmatrix
  from pyqmc.matrices.gms import Fort70

  f70 = Fort70(fort70file, nup=chkdata.nup, ndn=chkdata.ndn, verbose=verbose)

  PsiT = slaterdet.MultiDet()
  PsiT.dets = f70.psiT_det
  PsiT.set_ampl(f70.ampl)
  return PsiT



# The following routines are my bread-and-butter analysis tools, but they are not
# necessarily most general as the rest of the library.


def dets_calc_ovlp_with_PsiT1(wlkrs, psiT_det=None, psiT_kwd='PsiT'):
  """Calculate the overlap between a multideterminant wfn and a
  single-det wave function.
  Records the individual overlap on each determinant object for subsequent
  analysis."""
  from pyqmc.matrices.slaterdet import up_ovlp, dn_ovlp
  up_o_attr = 'up_%s_ovlp' % psiT_kwd
  dn_o_attr = 'dn_%s_ovlp' % psiT_kwd
  o_attr = '%s_ovlp' % psiT_kwd
  tot_ovlp = 0
  for D in wlkrs.dets:
    up_o, dn_o = up_ovlp(psiT_det, D), dn_ovlp(psiT_det, D)
    setattr(D, up_o_attr, up_o)
    setattr(D, dn_o_attr, dn_o)
    setattr(D, o_attr, up_o * dn_o)
    tot_ovlp += up_o * dn_o * D.ampl
  return tot_ovlp


def dets_calc_ovlp_with_PsiT_multidet(wlkrs, PsiT_det=None, psiT_kwd='PsiT'):
  """Calculate the overlap between a multideterminant wfn and a
  multidet wave function.
  Records the individual overlap on each determinant object for subsequent
  analysis."""
  from pyqmc.matrices.slaterdet import Det, Multidet, mdet_ovlp
  o_attr = '%s_ovlp' % psiT_kwd
  tot_ovlp = 0

  # Creates a dummy multidet object to work with mdet_ovlp below.
  D1 = Det()
  D1.ampl = 1.0
  MD = MultiDet()
  MD.dets = [D1]
  for D in wlkrs.dets:
    D1.up, D1.dn = D.up, D.dn
    ovlp_PsiT_phi = mdet_ovlp(PsiT, MD)
    setattr(D, o_attr, ovlp_PsiT_phi)
    tot_ovlp += ovlp_PsiT_phi * D.ampl
  return tot_ovlp


def dets_dump_stats(wlkrs, out=sys.stdout, psiT_kwd='PsiT'):
  """Verbosely reports the statistics of the determinants (weight,
  local energies, etc, and wave function overlap).
  """
  # In the development script elsewhere (Check_walkers.py),
  # it was called 'report_status_dets'.
  # But the layout has changed considerably.
  from wpylib.iofmt.text_output import text_output
  from wpylib.text_tools import str_fmt_heading
  fmt = "%4d  %17.9g  %17.9g  | %17.9g %17.9g | %17.9g %17.9gj | %17.9g %17.9gj | %17.9g %17.9gj || %17.9g %17.9gj | %17.9g %17.9gj | %17.9g %17.9gj |\n"
  fmt_heading = str_fmt_heading(fmt)
  #xkwd = lambda X : X.replace('%', psiT_kwd)
  heading = fmt_heading % ('no',
                           'wtwlkr', 'phasefac', 'Elocal', '(imag)',
                           'impfn', '(imag)',
                           'det_norm', '(imag)',
                           'ampl', '(imag)',
                           '%s_ovlp' % psiT_kwd, '(imag)',
                           'up_%s_ovlp' % psiT_kwd, '(imag)',
                           'dn_%s_ovlp' % psiT_kwd, '(imag)',
                          )
  Print = text_output(out)
  Print(heading)
  for (i, D) in enumerate(wlkrs.dets):
    impfn = getattr(D, 'impfn', 1e-99)
    det_norm = getattr(D, 'det_norm', 1e-99)
    ampl = getattr(D, 'ampl', 1e-99)
    o = getattr(D, '%s_ovlp' % psiT_kwd, 1e-99)
    up_o = getattr(D, 'up_%s_ovlp' % psiT_kwd, 1e-99)
    dn_o = getattr(D, 'dn_%s_ovlp' % psiT_kwd, 1e-99)

    Print(fmt \
          % (i+1,
             D.wtwlkr,
             D.phasefac,
             D.El.real, D.El.imag,
             impfn.real, impfn.imag,
             det_norm.real, det_norm.imag,
             ampl.real, ampl.imag,
             o.real, o.imag,
             up_o.real, dn_o.imag,
             dn_o.real, dn_o.imag,
            )
         )

