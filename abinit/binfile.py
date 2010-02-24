# $Id: binfile.py,v 1.3 2010-02-24 14:36:03 wirawan Exp $
#
# pyqmc.abinit.binfile
#
# Created: 20100208
# Wirawan Purwanto

import numpy
import sys

from wpylib.sugar import ifelse
from wpylib.iofmt.fortbin import fortran_bin_file

class abinit_header(object):
  pass

class abinit_psp_rec(object):
  pass

class abinit_bin_file(fortran_bin_file):
  supported_header = (42, 44, 53, 56, 57)
  int_ = numpy.int32
  double = numpy.float64

  #def open(self, filename, mode="r"):
  #  fortran_bin_file.open(self, filename, mode)

  def read_header(self):
    if self.debug:
      verbose = self.debug
      def dbg(lvl, msg):
        if verbose >= lvl: sys.stderr.write(msg)
    else:
      dbg = lambda lvl, msg : None

    self.header = abinit_header()
    hdr = self.header
    Int = self.int_
    Dbl = self.double
    self.read(('codvsn','S6'), ('headform',Int), ('fform',Int), dest=self)
    dbg(1, "codvsn=%s, headform=%d, fform=%d\n" \
        % (self.codvsn, self.headform, self.fform))

    dbg(1, "Reading header struct 1: bantot, date, ...\n")
    flds = self.get_header_structs(kind=1)
    self.read(dest=self.header, *flds)

    dbg(1, "Reading header struct 2: istwfk, nband...\n")
    flds2 = self.get_header_structs(kind=2)
    self.read(dest=self.header, *flds2)

    self.header.psp = []
    flds3 = self.get_header_structs(kind=3)
    for ipsp in xrange(hdr.npsp):
      newpsp = abinit_psp_rec()
      self.header.psp.append(newpsp)
      dbg(1, "Reading psp struct #%d\n" % ipsp)
      self.read(dest=newpsp, *flds3)

    dbg(1, "Reading header struct 4: residm, xred, ...\n")
    flds4 = self.get_header_structs(kind=4)
    self.read(dest=self.header, *flds4)

    # (in case of usepaw==1, there are some additional records)
    if getattr(hdr, "usepaw", None):
      raise NotImplementedError, \
        "usepaw!=0 has another header block; its reading is still unimplemented."

    self.header_fields = [ f[0] for f in flds + flds2 + flds4 ]
    self.header_fields_psp = [ f[0] for f in flds3 ]


  def write_header(self):
    """Puts a header back to disk with the data in the current object."""
    if self.debug:
      verbose = self.debug
      def dbg(lvl, msg):
        if lvl > verbose: sys.stderr.write(msg)
    else:
      dbg = lambda lvl, msg : None

    hdr = self.header
    Int = self.int_
    Dbl = self.double
    flds0 = [('codvsn','S6'), ('headform',Int), ('fform',Int)]
    self.write_fields(self, *flds0)

    dbg(1, "Writing header struct 1: bantot, date, ...\n")
    flds = self.get_header_structs(kind=1)
    self.write_fields(self.header, *flds)

    dbg(1, "Writing header struct 2: istwfk, nband...\n")
    flds2 = self.get_header_structs(kind=2)
    self.write_fields(self.header, *flds2)

    self.psp = []
    flds3 = self.get_header_structs(kind=3)
    for ipsp in xrange(hdr.npsp):
      newpsp = self.header.psp[ipsp]
      dbg(1, "Writing psp struct #%d\n" % ipsp)
      self.write_fields(newpsp, *flds3)

    dbg(1, "Writing header struct 4: residm, xred, ...\n")
    flds4 = self.get_header_structs(kind=4)
    self.write_fields(self.header, *flds4)

    # (in case of usepaw==1, there are some additional records)
    if getattr(hdr, "usepaw", None):
      raise NotImplementedError, \
        "usepaw!=0 has another header block; its writing is still unimplemented."


  def print_header(self, out=sys.stdout):
    out.write("codvsn=%s, headform=%d, fform=%d\n" \
              % (self.codvsn, self.headform, self.fform))
    def dump(varname):
      val = getattr(self.header, varname)
      if type(val) == self.int_:
        out.write("%s=%d\n" % (varname, val))
      elif type(val) == self.double:
        out.write("%s=%.15g\n" % (varname, val))
      else:
        out.write("%s=\n%s\n" % (varname, str(val)))
    #out.write("bantot=%d\n" % self.bantot)

    for f in self.header_fields:
      dump(f)
    #read(fdWFK) codvsn,headform,fform


  def load_density(self):
    """Loads the density data after the header.
    This is for the 'O_DEN' output file.
    Only supports real densities right now. (FIXME)

    """
    hdr = self.header
    fftsize = (hdr.ngfft[0],hdr.ngfft[1],hdr.ngfft[2])
    self.read(('rhor', self.double, fftsize), dest=self)

    if hdr.nspden == 1:
      pass
    elif hdr.nspden == 2:
      self.read(('rhor_up', self.double, fftsize), dest=self)
      self.rhor_dn = self.rhor - self.rhor_up
    else:
      raise NotImplementedError, \
        "Density reading is not implemented completely for nspden==%d" % hdr.nspden

    return self.rhor

  def get_header_structs(self, kind):
    """Obtains various ABINIT header structures; one at a time.
    Caveat: Sizes in some header structures depend on the preceding structure."""

    Int = self.int_
    Dbl = self.double
    hdr = self.header
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

    if kind == 1: return [
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

    if kind == 2: return [
      ('istwfk', Int, hdr.nkpt),
      ('nband', Int, (hdr.nkpt,hdr.nsppol)),
      ('npwarr', Int, (hdr.nkpt)),
    ] \
    + ifelse(self.headform <= 44, [('so_typat', Int, (hdr.ntypat))], []) \
    + ifelse(self.headform >= 53, [('so_psp', Int, (hdr.npsp))], []) \
    + [
      ('symafm', Int, (hdr.nsym)),
      ('symrel', Int, (3,3,hdr.nsym)),
      ('typat', Int, (hdr.natom)),
      ('kpt', Dbl, (3,hdr.nkpt)),
      ('occ', Dbl, (hdr.bantot)),
      ('tnons', Dbl, (3,hdr.nsym)),
      ('znucltypat', Dbl, (hdr.ntypat)),
    ] \
    + ifelse(self.headform >= 53, [('wtk', Dbl, (hdr.nkpt))], []) \
    + [
    ]

    """
  do ipsp=1,npsp
  ! (npsp lines, 1 for each pseudopotential; npsp=ntypat, except if
  ! alchemical pseudo-atoms)
<= 4.2:
      read (fdWFK) title,znuclpsp,zionpsp,pspso,pspdat,pspcod,pspxc

>= 4.4:
      read (fdWFK) title,znuclpsp,zionpsp,pspso,pspdat,pspcod,pspxc,lmn_size

    """
    if kind == 3: return [ # psp header
      ('title', 'S132'),
      ('znuclpsp', Dbl),
      ('zionpsp', Dbl),
      ('pspso', Int),
      ('pspdat', Int),
      ('pspcod', Int),
      ('pspxc', Int),
    ] \
    + ifelse(self.headform >= 44, [('lmn_size', Int)], []) \
    + [
    ]

    """
(in case of usepaw==0, final record: residm, coordinates, total energy,
Fermi energy)

      write(unit=unit) residm,xred(1:3,1:natom),etotal,fermie

    """
    if kind == 4: return [
      ('residm', Dbl),
      ('xred', Dbl, (3, hdr.natom)),
      ('etotal', Dbl),
      ('fermie', Dbl),
    ]

    """
(in case of usepaw==1, there are some additional records)

    """
    raise ValueError, "Invalid header structure kind = %d" % kind
    #return (flds, flds2, flds3, flds4)


class abinit_density_file(abinit_bin_file):
  pass

