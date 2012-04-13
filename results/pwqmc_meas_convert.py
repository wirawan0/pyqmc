#!/usr/bin/python
# $Id: pwqmc_meas_convert.py,v 1.1 2010-09-30 01:46:43 wirawan Exp $
#
# pwqmc_meas_convert.py
# Tools to convert PWQMC-77 measurement dump (pwaf-*.ene) files
# to the standard HDF5 format.
#
# Wirawan Purwanto
# Created: 20100513
#

import os
import time

import h5py # HDF5 Python API
import numpy

import pyqmc.results.pwqmc_info as pwinfo
import pyqmc.results.pwqmc_meas as pwmeas

from pyqmc.results.pwqmc_info import \
  kptstr, pwqmc_info

from wpylib.timer import timer
from wpylib.sugar import Parameters


class convert_ene2hdf5(object):
  """Main converter from *.ene files to a HDF5-formatted measurement file.
  Most likely you will need to derive this class to make a converter suitable
  to your own need."""
  def __init__(self, **opts):
    self.opts = Parameters()
    self.opts.TMPDIR = sh.getenv("PYQMC_TMPDIR", "TMPDIR", "TMP", default="/tmp") \
        + ("/pwqmc_convert.%d" % (os.getpid()))
    self.opts.cleanup_tmpdir = True

  """
  def get_opt_proc(self, opts):
    def_opts = self.opts
    def get(opt, default):
      if opt in opts:
        return opts[opt]
      elif kw in
        getattr(opts, kw)
  """

  def convert(self, src=None, info=None, output=None,
              _defaults=dict(src="PWAF-meas.tar.lzma", info="INFO",
                             output="measurements.h5", hints={},
                             E_prefactor=1.0,
                             debug=1),
              **_opts_):
    #def convert(src="PWAF-meas.tar.lzma", info="INFO",
    #            output="measurements.h5",
    #            db_hints={},
    #            debug=1):
    """Converts a set of measurement data (*.ene) to a standard meas_hdf5's
    HDF5 database.
    Also adds some useful metadata from the INFO file.
    Required parameters to be properly set:
    * src = the tarball, or the glob of the measurement files (*.ene).
    * info = the INFO file.

    Additional methods, if defined in the derived object, will be called:
    * convert_preamble_steps_(hint, info, opts)
    * convert_postamble_steps_(hdf5_raw_group, info, opts)
    """
    # FIXME: use self-introspection to reduce kitchen-sink params here:
    #p = Parameters(locals(), _opts_, _opts_.get('opts'), self.opts, _defaults)
    p = opts._create_(_defaults)
    info_file = p.info
    src = p.src
    orig_dir = os.getcwd()

    tm1 = time.clock()
    if not isinstance(info_file, pwqmc_info):
      info = pwqmc_info(info_file)
    else:
      info = info_file
      info_file = info['info_file']

    kpt_str = kptstr(info['kpt'])
    # do these fetches here just in case they fail due to my mistake/negligence
    #kpt_data1 = ALL_KPTS_DATA[cellstr][volstr][kpt_str]
    #E_GGA = ALL_KPTS_DATA[cellstr][volstr]['kgrid']['E_GGA']
    if 'Etrial_noconst' not in info:
      raise RuntimeError, \
        "Trial energy is not found in the INFO file: %s" % (info_file)

    use_tmpdir = False
    if isinstance(src, basestring):
      if src.endswith(".tar.lzma"):
        os.chdir(p.TMPDIR)
        sh.system("rm -rf *" % (p.TMPDIR,)) # always start clean
        os.chdir(orig_dir)
        os.run("tar", ("-C", p.TMPDIR, "--use-compress-program=lzma", "-xf", src))
        files = TMPDIR + "/*.ene"
        use_tmpdir = True
      elif src.endswith(".tar.bz2"):
        os.chdir(p.TMPDIR)
        sh.system("rm -rf *" % (p.TMPDIR,)) # always start clean
        os.chdir(orig_dir)
        os.run("tar", ("-C", p.TMPDIR, "-j", "-xf", src))
        files = TMPDIR + "/*.ene"
        use_tmpdir = True
      elif src.endswith(".ene"): # or src.endswith(".meas"):
        # *.meas is OLD. Don't use anymore.
        files = src
      else:
        raise ValueError, "Don't know how to handle src = %s" % (src)
    else:
      raise ValueError, "Don't know how to handle src %s = %s" % (str(type(src)), str(src))

    hints = hints.copy()
    hints.update({
      'nwlkavg': info['nwlk'],
      'nwlkmax': info['nwlkmax'],
      #'default_raw_chunks': [1, info['nwlkmax']],
      #'value_processor': valpx,
    })

    if p.E_prefactor != 1.0 and 'value_processor' not in hints:
      def valpx(data, meta, *junk1, **junk2):
        """Renormalize the energy values (real and imaginary!) to 4-atom
        cell value."""
        data['E_l'] *= E_prefactor
      hints['value_processor'] = valpx

    if 'convert_preamble_steps_' in dir(self):
      # This is useful for e.g. adding default_raw_chunks, defining
      # value_processor, etc.
      self.convert_preamble_steps_(hints=hints, info=info, opts=p)
      # Examples:
      # in MnO 2x2x2:
      #    natoms = cell_info[cellstr]['natoms']
      #    E_prefactor = 4.0 / natoms
      #    def valpx(data, meta, *junk1, **junk2):
      #      """Renormalize the energy values (real and imaginary!) to 4-atom
      #      cell value."""
      #      data['E_l'] *= E_prefactor
      # NOTE: special to MnO runs,
      # We use chunksize = [1,nwlkmax] because we know the popsizes are
      # hovering near nwlkmax anyway.

    db = pwmeas.convert_meas_to_hdf5_v2(output,
           files=files,
           betablk=info["betablk"], deltau=info["deltau"],
           H0=info["H0"],
           debug=debug,
           **hints
         )
    db.flush()
    if debug > 0:
      self.last_db = db

    # The last opened raw group is given as db.raw
    # --this is also the data group which contain the dataset from the
    # conversion just done.
    job0 = db.raw.job()
    # Add some useful attributes:
    job0.attrs['System'] = info['System']
    job0.attrs['H0'] = info['H0']
    job0.attrs['deltau'] = info['deltau']
    job0.attrs['Evar'] = info['Evar'] * E_prefactor
    ET = (info['Etrial_noconst'] + info['H0']) * E_prefactor
    ET_delta = 2 / math.sqrt(info['deltau']) * E_prefactor
    job0.attrs['Etrial'] = ET # for El_bounds
    job0.attrs['Ebounds'] = (ET - ET_delta, ET + ET_delta)
    job0.attrs['units'] = 'Ry'  # temporarily storing in Rydberg unit
    job0.attrs['kpt'] = info['kpt']
    job0.attrs['kptstr'] = "k" + kpt_str # compatibility with old format
    job0.attrs['wsvol'] = info['vol']
    #vol_code = "%.2f" % (info['vol'] * bohr__angstrom**3 / 2)
    #job0.attrs['GeomCode'] = "vol" + vol_code

    # Additional metadata can be specified afterward
    if 'System' in p:
      System = p.System
    elif 'System' in info:
      System = info['System']

    if 'convert_postamble_steps_' in dir(self):
      # This is useful for e.g. adding more metadata
      self.convert_postamble_steps_(raw_group=job0, info=info, opts=p)

    tm2 = time.clock()
    print "%s : kpt_str = %s, ET = %g, E_GGA_kpt = %g, E_GGA = %g; total time = %d" % \
      (output, job0.attrs['kptstr'],
       job0.attrs['Etrial'], job0.attrs['E_GGA_kpt'], job0.attrs['E_GGA'],
       tm2 - tm1
      )
    db.flush()
    return db

