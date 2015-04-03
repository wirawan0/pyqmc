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
  """
  from numpy import zeros, eye, asmatrix
  from pyqmc.matrices.gms import Fort70

  f70 = Fort70(fort70file, nup=chkdata.nup, ndn=chkdata.ndn, verbose=verbose)

  PsiT = slaterdet.MultiDet()
  PsiT.dets = f70.psiT_det
  PsiT.set_ampl(f70.ampl)
  return PsiT


