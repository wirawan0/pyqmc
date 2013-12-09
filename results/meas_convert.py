#!/usr/bin/python
# $Id: pwqmc_meas_convert.py,v 1.1 2010-09-30 01:46:43 wirawan Exp $
#
# pyqmc.results.meas_convert.py
# Tools to convert PWQMC-77/GAFQMC measurement dump (pwaf-*.ene) files
# to the standard HDF5 format.
#
# Wirawan Purwanto
# Created: 20100513
#

import math
import os
import time

import h5py # HDF5 Python API
import numpy

from copy import copy

import pyqmc.results.pwqmc_info as pwinfo
import pyqmc.results.pwqmc_meas as pwmeas
import pyqmc.results.gafqmc_info as ginfo
import pyqmc.results.gafqmc_meas as gmeas

from pyqmc.results.pwqmc_info import \
  kptstr, pwqmc_info

from pyqmc.results.gafqmc_info import \
  gafqmc_info

from wpylib.timer import timer
from wpylib.text_tools import str_trunc_begin
from wpylib.params.params_flat import Parameters
import wpylib.shell_tools as sh


def store_pwqmc_info_metadata(info, hdf5_job_rec):
  """Stores additional, basis/method-specific metadata on
  the HDF5 job record.

  This is for PWQMC runs."""
  job0 = hdf5_job_rec

  kpt_str = kptstr(info['kpt'])
  job0.attrs['kpt'] = info['kpt']
  job0.attrs['kptstr'] = "k" + kpt_str # compatibility with old format
  job0.attrs['wsvol'] = info['vol']


def store_gafqmc_info_metadata(info, hdf5_job_rec):
  """Stores additional, basis/method-specific metadata on
  the HDF5 job record.

  This is for GAFQMC runs."""
  pass


def dict_nested_copy_two_levels(d):
  """Do dict nested copy for two outermost levels;
  beyond that no copy is made."""
  c = {}
  for k in d:
    c[k] = d[k].copy()
  return c

"""
Default parameters for convert:

db_hints: a dict, containing default parameters to be passed on to
the actual converter machinery in pwqmc_meas module.

"""

class convert_ene2hdf5(object):
  """Main converter from *.ene files to a HDF5-formatted measurement file.
  Most likely you will need to derive this class to make a converter suitable
  to your own need."""

  # Change these if needed:
  meas_module = pwmeas
  info_module = pwinfo
  info_class = pwinfo.pwqmc_info

  Default_desc = {
    # Metadata descriptor stored in the HDF5 file.
    # System = descriptive text of the physical system
    # unit = units of the energy measurement
    # Note: This object will be nested-copied by two levels upon the
    # instantiation of the object.
    pwinfo.pwqmc_info: dict(
      System="PWQMC calculation",
      unit="Ry",
      extra_meta_copy=store_pwqmc_info_metadata,
    ), # pwqmc_info
    ginfo.gafqmc_info: dict(
      System="GAFQMC calculation",
      unit="Ha",
      extra_meta_copy=store_gafqmc_info_metadata,
    ), # gafqmc_info
  }

  Default_params = {
    # User-overridable parameters for function calls below.
    # Note: This object will be nested-copied by two levels upon the
    # instantiation of the object.
    pwinfo.pwqmc_info: dict(
      convert_defaults=dict(
        src="PWAF-meas.tar.lzma",
        info="INFO",
        output="measurements.h5",
        db_hints=dict(keep_phasefac='auto',keep_El_imag='auto',),
        E_prefactor=1.0,
        debug=1
      ),
    ), # pwqmc_info
    ginfo.gafqmc_info: dict(
      convert_defaults=dict(
        src="GAFQMCF-meas.tar.lzma",
        info="INFO",
        output="measurements.h5",
        db_hints=dict(keep_phasefac='auto',keep_El_imag='auto',),
        E_prefactor=1.0,
        debug=1
      ),
    ), # gafqmc_info
  }

  def __init__(self, **opts):
    self.opts = Parameters()
    self.opts.TMPDIR = sh.getenv("PYQMC_TMPDIR", "TMPDIR", "TMP", default="/tmp") \
        + ("/pwqmc_convert.%d" % (os.getpid()))
    self.opts.cleanup_tmpdir = True
    self.Default_params = dict_nested_copy_two_levels(self.Default_params)
    self.Default_desc = dict_nested_copy_two_levels(self.Default_desc)

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
              **_opts_):
    #def convert(src="PWAF-meas.tar.lzma", info="INFO",
    #            output="measurements.h5",
    #            db_hints={},
    #            debug=1):
    """Converts a set of measurement data (*.ene) to a standard meas_hdf5's
    HDF5 database.
    Also adds some useful metadata from the INFO file.
    Required parameters to be properly set:
    * src = the filename of the tarball, or the glob of
      the measurement files (*.ene).
    * info = the INFO file.
    * output = the target HDF5 data file.

    Additional methods, if defined in the derived object, will be called:
    * convert_preamble_steps_(hint, info, opts)
    * convert_postamble_steps_(hdf5_raw_group, info, opts)
    """
    # FIXME: use self-introspection to reduce kitchen-sink params here:
    #p = Parameters(locals(), _opts_, _opts_.get('opts'), self.opts, _defaults)
    # The function defaults are now provided in the
    # Default_params...['convert_defaults'] field.
    p = self.opts._create_(self.Default_params[self.info_class]['convert_defaults'])
    if info == None:
      info_file = p.info
    else:
      info_file = info
    if src == None: src = p.src
    if output == None: output = p.output
    orig_dir = os.getcwd()

    tm1 = time.clock()
    if not isinstance(info_file, pwqmc_info):
      info = self.info_class(info_file)
    else:
      info = info_file
      info_file = info['info_file']

    # Deduce the datatype of the info structure:
    info_dtype = None
    for klass in self.Default_desc.keys():
      if isinstance(info, klass):
        info_dtype = klass
        break
    if info_dtype == None:
      raise RuntimeError, \
        "Cannot deduce the datatype of the info structure: %s" % type(info)

    Default_desc = self.Default_desc[info_dtype]
    if not issubclass(info_dtype, self.info_class):
      from warnings import warn
      warn("INFO class (%s) does not match the converter's info_class!" \
           % (info_dtype, self.info_class),
           UserWarning)

    # do these fetches here just in case they fail due to my mistake/negligence
    #kpt_data1 = ALL_KPTS_DATA[cellstr][volstr][kpt_str]
    #E_GGA = ALL_KPTS_DATA[cellstr][volstr]['kgrid']['E_GGA']
    if 'Etrial_noconst' not in info:
      raise PyqmcDataError, \
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

    try:
      db_hints = dict(p.db_hints)
    except:
      db_hints = {}

    db_hints.update({
      'nwlkavg': info['nwlk'],
      'nwlkmax': info['nwlkmax'],
      #'default_raw_chunks': [1, info['nwlkmax']],
      #'value_processor': valpx,
    })
    is_free_proj = info['constraint'] in ('none',)
    if db_hints.get('keep_El_imag') == 'auto':
      db_hints['keep_El_imag'] = is_free_proj
    if db_hints.get('keep_phasefac') == 'auto':
      db_hints['keep_phasefac'] = is_free_proj

    if p.E_prefactor != 1.0 and 'value_processor' not in db_hints:
      def valpx(data, meta, *junk1, **junk2):
        """Rescales the energy values (real and imaginary!)."""
        data['E_l'] *= p.E_prefactor
      db_hints['value_processor'] = valpx

    if 'convert_preamble_steps_' in dir(self):
      # This is useful for e.g. adding default_raw_chunks, defining
      # value_processor, etc.
      self.convert_preamble_steps_(hints=db_hints, info=info, opts=p)
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

    db = self.meas_module.convert_meas_to_hdf5_v2(output,
           files=files,
           betablk=info["betablk"], deltau=info["deltau"],
           H0=info["H0"],
           debug=p.debug,
           **db_hints
         )
    db.flush()
    if p.debug > 0:
      self.last_db = db

    # The last opened raw group is given as db.raw
    # --this is also the data group which contain the dataset from the
    # conversion just done.
    job0 = db.raw.job()
    # Add some useful attributes:
    if 'System' in p:
      System = p.System
    elif 'System' in info:
      System = info['System']
    elif 'System' in Default_desc:
      System = Default_desc['System']
    else:
      System = '(Unknown calculation)'
    job0.attrs['System'] = System
    job0.attrs['H0'] = info['H0']
    job0.attrs['deltau'] = info['deltau']
    job0.attrs['Evar'] = info['Evar'] * p.E_prefactor
    ET = (info['Etrial_noconst'] + info['H0']) * p.E_prefactor
    ET_delta = 2 / math.sqrt(info['deltau']) * p.E_prefactor
    job0.attrs['Etrial'] = ET * p.E_prefactor # for El_bounds
    job0.attrs['Ebounds'] = ((ET - ET_delta) * p.E_prefactor, (ET + ET_delta) * p.E_prefactor)
    job0.attrs['units'] = Default_desc['unit']
    extra_meta_copy = Default_desc.get('extra_meta_copy', None)
    if extra_meta_copy != None:
      extra_meta_copy(info, job0)

    if 'convert_postamble_steps_' in dir(self):
      # This is useful for e.g. adding more metadata
      self.convert_postamble_steps_(raw_group=job0, info=info, opts=p)

    tm2 = time.clock()
    print "%s : ET = %g; total time = %d; %s" % \
      (output,
       job0.attrs['Etrial'],
       tm2 - tm1,
       str_trunc_begin(System, 64),
      )

    db.flush()
    return db

  @property
  def convert_defaults(self):
    return self.Default_params[self.info_class]['convert_defaults']


class gafqmc_convert_ene2hdf5(convert_ene2hdf5):
  """An energy measurement data format converter for PWQMC calculations.
  """
  def __init__(self, **opts):
    """Creates a energy measurement data format converter.
    For GAFQMC calculations.
    """
    self.meas_module = gmeas
    self.info_module = ginfo
    self.info_class = ginfo.gafqmc_info
    super(gafqmc_convert_ene2hdf5, self).__init__(**opts)


class pwqmc_convert_ene2hdf5(convert_ene2hdf5):
  """An energy measurement data format converter for PWQMC calculations.
  """
  pass

