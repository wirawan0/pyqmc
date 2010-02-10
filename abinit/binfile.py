# $Id: binfile.py,v 1.1 2010-02-10 23:11:55 wirawan Exp $
#
# pyqmc.abinit.binfile
#
# Created: 20100208
# Wirawan Purwanto

import numpy
import sys

from wpylib.sugar import ifelse
from wpylib.iofmt.fortbin import fortran_bin_file

class abinit_bin_file(fortran_bin_file):
  supported_header = (42, 44, 53, 56, 57)
  int = numpy.int32
  double = numpy.float64
  def open(self, filename):
    self.F = open(filename, "rb")
  def read_header(self, verbose=100):
    if verbose:
      def dbg(lvl, msg):
        if lvl > verbose: sys.stderr.write(msg)
    else:
      dbg = lambda lvl, msg : None

    Int = self.int
    Dbl = self.double
    self.read(('codvsn','S6'), ('headform',Int), ('fform',Int), out=self)
    dbg(1, "codvsn=%s, headform=%d, fform=%d\n" \
        % (self.codvsn, self.headform, self.fform))

    """From ABINIT manuals:

4.2:

    read(fdWFK) bantot,date,intxc,ixc,natom,ngfft(1:3), &
     & nkpt,nspden,nspinor,nsppol,nsym,npsp,ntypat,occopt,pertcase, &
     & ecut,ecutsm,ecut_eff,qptn(1:3),rprimd(1:3,1:3),stmbias,tphysel,tsmear


4.4, 5.3, 5.6:

    read(fdWFK) bantot,date,intxc,ixc,natom,ngfft(1:3), &
     & nkpt,nspden,nspinor,nsppol,nsym,npsp,ntypat,occopt,pertcase, &
     & usepaw, &
     & ecut,ecutdg,ecutsm,ecut_eff,qptn(1:3), &
     & rprimd(1:3,1:3),stmbias,tphysel,tsmear

5.7:

    read(fdWFK) bantot,date,intxc,ixc,natom,ngfft(1:3), &
     & nkpt,nspden,nspinor,nsppol,nsym,npsp,ntypat,occopt,pertcase, &
     & usepaw, &
     & ecut,ecutdg,ecutsm,ecut_eff,qptn(1:3), &
     & rprimd(1:3,1:3),stmbias,tphysel,tsmear,usewvl

    """

    flds = [
      ('bantot', Int),
      ('date', Int),
      ('intxc', Int),
      ('ixc', Int),
      ('natom', Int),
      ('ngfft', Int, 3),
      ('nkpt', Int),
      ('nspden', Int),
      ('nspinor', Int),
      ('nsppol', Int),
      ('nsym', Int),
      ('npsp', Int),
      ('ntypat', Int),
      ('occopt', Int),
      ('pertcase', Int),
    ] \
    + ifelse(self.headform >= 44, [('usepaw', Int)], []) \
    + [
      ('ecut', Dbl),
    ] \
    + ifelse(self.headform >= 44, [('ecutdg', Dbl)], []) \
    + [
      ('ecutsm', Dbl),
      ('ecut_eff', Dbl),
      ('qptn', Dbl, 3),
      ('rprimd', Dbl, (3,3)),
      ('stmbias', Dbl),
      ('tphysel', Dbl),
      ('tsmear', Dbl),
    ] \
    + ifelse(self.headform >= 57, [('usewvl', Int)], []) \
    + [
    ]
    self.read(*flds, out=self)

    """
<= 4.4:

    read(fdWFK) istwfk(1:nkpt), nband(1:nkpt*nsppol), &
     & npwarr(1:nkpt), so_typat(1:ntypat), &
     & symafm(1:nsym), &
     & symrel(1:3,1:3,1:nsym), &
     & typat(1:natom), kpt(1:3,1:nkpt), &
     & occ(1:bantot), tnons(1:3,1:nsym), &
     & znucltypat(1:ntypat)

>= 5.3:

    read(fdWFK) istwfk(1:nkpt), nband(1:nkpt*nsppol), &
     & npwarr(1:nkpt), so_psp(1:npsp), &
     & symafm(1:nsym), &
     & symrel(1:3,1:3,1:nsym), &
     & typat(1:natom), kpt(1:3,1:nkpt), &
     & occ(1:bantot), tnons(1:3,1:nsym), &
     & znucltypat(1:ntypat), &
     & wtk(1:nkpt)
    """

    flds2 = [
      ('istwfk', Int, self.nkpt),
      ('nband', Int, (self.nkpt,self.nsppol)),
      ('npwarr', Int, (self.nkpt)),
    ] \
    + ifelse(self.headform <= 44, [('so_typat', Int, (self.ntypat))], []) \
    + ifelse(self.headform >= 53, [('so_psp', Int, (self.npsp))], []) \
    + [
      ('symafm', Int, (self.nsym)),
      ('symrel', Int, (3,3,self.nsym)),
      ('typat', Int, (self.natom)),
      ('kpt', Dbl, (3,self.nkpt)),
      ('occ', Dbl, (self.bantot)),
      ('tnons', Dbl, (3,self.nsym)),
      ('znucltypat', Dbl, (self.ntypat)),
    ] \
    + ifelse(self.headform >= 53, [('wtk', Dbl, (self.nkpt))], []) \
    + [
    ]
    print repr(flds2)
    self.read(*flds2, out=self)

    self.header_fields = [ f[0] for f in flds + flds2 ]


  def print_header(self, out=sys.stdout):
    out.write("codvsn=%s, headform=%d, fform=%d\n" \
              % (self.codvsn, self.headform, self.fform))
    def dump(varname):
      val = getattr(self, varname)
      if type(val) == self.int:
        out.write("%s=%d\n" % (varname, val))
      elif type(val) == self.double:
        out.write("%s=%.15g\n" % (varname, val))
      else:
        out.write("%s=\n%s\n" % (varname, str(val)))
    #out.write("bantot=%d\n" % self.bantot)

    for f in self.header_fields:
      dump(f)
    #read(fdWFK) codvsn,headform,fform



class abinit_density_file(abinit_bin_file):
  pass

