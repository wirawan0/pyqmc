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

  # By default assume the walkers reside in a subdir called "wlk/"
  def read(self, src, verbose=0, output=sys.stdout):
    """Opens the `src' file and loads the walker data."""
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
    F = fortran_bin_file(src)
    F.read(('code_version', Float), dest=rec)
    F.read(('date', Str(10)), dest=rec)
    F.read(('time', Str(10)), dest=rec)
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
    F.read(('fmt_version', Float), dest=rec) # checkpoint format version
    F.read(('nh', Int), ('anorm', Float), ('etrial', Float),
           ('istpacc', Int), dest=rec)
    F.read(('timeshift', Float), dest=rec)
    F.read(('nwlkr_proc', Int), dest=rec)
    if verbose >= 10:
      w(("  lran        = %s\n" % rec.lran) + \
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
      w("# Walker data: ampl, phasefac, re(El), im(El)\n")

    pop = MultiDet()
    pop.dets = []
    rec.wlkrs = pop
    for iwlk in xrange(rec.nwlkr_proc):
      wlk = Det()
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


