#!/usr/bin/python
# $Id: gafqmc_walkers.py,v 1.1 2011-09-01 20:12:49 wirawan Exp $
#
# pyqmc.data.gafqmc_walkers module
#
# Tools for handling GAFQMC walker files
#
# Wirawan Purwanto
# Created: 20090320
#

import numpy
import sys

from wpylib.params.params_struct import Struct as struct

from wpylib.iofmt.fortbin import fortran_bin_file
from wpylib.iofmt.text_output import text_output

from pyqmc.matrices.slaterdet import Det, MultiDet

class wlk_gafqmc(object):
  """Legacy GAFQMC-style walker files.
  Each MPI process has its own W_old_run_NNN file.
  We only support serial version of walker files.

  Internally, the walkers are always stored in "unrestricted" format
  (see pyqmc.matrices.slaterdet.Det).

  Potential problem with MPI's W_old_run_NNN files:
    * in the QMC code there is lran1..4 for myid > 0
      only if fast==1. This causes the format to be inconsistent,
      depending on run configuration.
  """
  Complex = fortran_bin_file.default_complex
  Float = fortran_bin_file.default_float
  Int = fortran_bin_file.default_int
  def init_metadata(self, nbasis, nup, ndn, udet):
    """Initializes some important metadata for the walker file.
    Unfortunately this set of data is not provided inside the walker
    file itself.
    """
    self.nbasis = nbasis
    self.nup = nup
    self.ndn = ndn
    self.udet = udet
    assert nup >= ndn

  def read(self, src, verbose=0, output=sys.stdout, rank=0, indep_rng=1):
    """Opens the `src' file and loads the walker data.
    Input parameters:
    - verbose = log verbosity level (default 0)
    - output = file-like stream for log output (default stdout)
    - rank = MPI rank of the process producing this checkpoint file
    - indep_rng = indicator whether an independent random number generator
      is used per each MPI process.
    """
    Complex = self.Complex
    Float = self.Float
    Int = self.Int
    Str = lambda length : numpy.dtype('S' + str(length))
    if verbose:
      w = text_output(output, flush=True)
    else:
      w = text_output(None)
    rec = struct()
    self.data = rec
    w("Reading checkpoint file %s\n" % (src,))
    F = fortran_bin_file(src)
    if rank == 0:
      F.read(('code_version', Float), dest=rec)
      F.read(('date', Str(10)), dest=rec)
      F.read(('time', Str(10)), dest=rec)
      if verbose >= 10:
        w(("GAFQMC walker file data\n" \
           "  code_version = %.14g\n" \
           "  date = %s\n" \
           "  time = %s\n") \
           % (rec.code_version, rec.date, rec.time))
      F.read(('lran', Int, 4), dest=rec)
      F.read(('nblktot', Int), dest=rec)
      F.read(('uptot', Float), dest=rec)
      F.read(('downtot', Float), dest=rec)
      F.read(('srun', Float), dest=rec)
      F.read(('s2run', Float), dest=rec)
    elif indep_rng:
      F.read(('lran', Int, 4), dest=rec)
    F.read(('fmt_version', Float), dest=rec) # checkpoint format version
    F.read(('nh', Int), ('anorm', Float), ('etrial', Float),
           ('istpacc', Int), dest=rec)
    F.read(('timeshift', Float), dest=rec)
    F.read(('nwlkr_proc', Int), dest=rec)
    if verbose >= 10:
      w(("  lran        = %s\n" % rec.lran))
      if rank == 0:
        w( \
          ("  nblktot     = %i\n" % rec.nblktot) + \
          ("  uptot       = %.14g\n" % rec.uptot) + \
          ("  downtot     = %.14g\n" % rec.downtot) + \
          ("  srun        = %.14g\n" % rec.srun) + \
          ("  s2run       = %.14g\n" % rec.s2run) + \
          ("  fmt_version = %.14g\n" % rec.fmt_version) + \
          ("  nh          = %i\n" % rec.nh) + \
          ("  anorm       = %.14g\n" % rec.anorm) + \
          ("  etrial      = %.14g\n" % rec.etrial) + \
          ("  istpacc     = %i\n" % rec.istpacc) + \
          ("  timeshift   = %i\n" % rec.timeshift) + \
          ("  nwlkr_proc  = %i\n" % rec.nwlkr_proc) + \
        "")

    if self.udet:
      nptot = self.nup + self.ndn
    else:
      nptot = self.nup

    if verbose >= 20:
      w("# Walker data: wtwlkr, phasefac, re(El), im(El)\n")

    pop = MultiDet()
    pop.dets = []
    rec.wlkrs = pop
    for iwlk in xrange(rec.nwlkr_proc):
      wlk = Det()
      wlk.proc_rank = rank
      pop.dets.append(wlk)
      F.read(('iw', Int), dest=wlk)
      F.read(('wtwlkr', Float), dest=wlk)
      #wlk.wtwlkr = wlk.ampl
      # Note: the correct ampl is wtwlkr / <psiT|wlk> .
      # This must be applied later.
      F.read(('phasefac', Float), dest=wlk)
      F.read(('El', Complex), dest=wlk)
      orbs = numpy.zeros((self.nbasis, nptot), dtype=Complex)
      for ip in xrange(nptot):
        orbs[:,ip] = F.read(('col', Complex, self.nbasis))['col']
      if self.udet:
        wlk.make(src_up=orbs[:,0:self.nup], src_dn=orbs[:,self.nup:nptot])
      else:
        wlk.make(src_up=orbs[:,0:self.nup], src_dn=orbs[:,0:self.ndn])
      if verbose >= 20:
        w("%5d %14.10f %14.10f %16.10f %14.10f\n" \
          % (wlk.iw, wlk.wtwlkr, wlk.phasefac, wlk.El.real, wlk.El.imag))
    w("  :: %d walkers read from file %s\n" % (len(pop), src))

    if rec.fmt_version >= 2.0:
      F.read(('iflg_chkpt_impfn', Int), dest=rec)
      if rec.iflg_chkpt_impfn > 0:
        F.read(('nwlkr_proc_1', Int), dest=rec)
        assert rec.nwlkr_proc_1 == rec.nwlkr_proc
        F.read(('impfn', Complex, (rec.nwlkr_proc_1,)), dest=rec)
        # todo: affix impfn to Det objects above
        # NOTE: if you use read_all() method instead, these values
        # would have been affixed there!

    return rec


  def read_all(self, procs, fname_pattern='wlk/gafqmc-%(rank)05d',
               verbose=0, output=sys.stdout, indep_rng=1, fname_args={}):
    """Reads all walker files into a single big result record.
    The `procs` argument can be one of the following:
    - an integer > 0, which is the number of MPI processes
      (for parallel run), indicating that all the walker files woud
      be read
    - a list/tuple/array of process ranks, from which associated
      checkpoint files we want to read in the walkers.

    By default, we assume that all the walker files reside in a subdir
    called "wlk/" .
    """
    from wpylib.sugar import is_iterable
    from pyqmc import PyqmcDataError
    from itertools import izip

    if verbose:
      w = text_output(output, flush=True)
    else:
      w = text_output(None)

    if not is_iterable(procs):
      ranks = xrange(procs)
    else:
      ranks = procs

    w("read_all: Reading from %d checkpoint files\n" % len(ranks))
    dest = None
    impfn_all = []
    has_impfn = None
    for rank in ranks:
      filename = fname_pattern % (dict(fname_args, rank=rank))
      chk = self.read(src=filename, verbose=verbose, output=w,
                      rank=rank, indep_rng=indep_rng)
      # don't leave the `data` field, or else it will confuse user later.
      # Use `data_all` instead!
      del self.data
      if dest is None:
        dest = chk
        self.data_all = chk
        dest.NUM_WARNINGS = 0
        dest.ranks = [ rank ]
        dest.lran_proc = { rank : dest.lran }
        del dest.lran
        dest.wlkrs_proc = { rank : dest.wlkrs }
        has_impfn = hasattr(chk, 'impfn')
      else:
        dest.ranks.append(rank)
        # Do some sanity checks and issues warning irregular stuff.
        def check_param(name, val, ref_val):
          if ref_val != val:
            w(" Warning: parameter `%s' is different from expected value (%s, ref: %s)\n" \
              % (name, val, ref_val))
            dest.NUM_WARNINGS = dest.NUM_WARNINGS + 1

        check_param('nh', chk.nh, dest.nh)
        check_param('anorm', chk.anorm, dest.anorm)
        check_param('etrial', chk.etrial, dest.etrial)
        check_param('istpacc', chk.istpacc, dest.istpacc)
        check_param('timeshift', chk.timeshift, dest.timeshift)

        dest.lran_proc[rank] = chk.lran
        dest.wlkrs_proc[rank] = chk.wlkrs
        dest.wlkrs.dets.extend(chk.wlkrs.dets)

        if has_impfn != hasattr(chk, 'impfn'):
          raise PyqmcDataError, \
            (("Inconsistent existence of impfn field across walker files " \
              " (currently on rank #%s; first-rank impfn status was %s)") \
             % rank)

      if has_impfn:
        impfn_all.append(chk.impfn)

    # Final brush-up:
    if has_impfn:
      dest.impfn = numpy.concatenate(impfn_all)
      # Affix impfn to D
      for (impfn, D) in izip(dest.impfn, dest.wlkrs):
        D.impfn = impfn

    return dest
