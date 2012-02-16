#!/usr/bin/python
# $Id: __init__.py,v 1.3 2011-09-12 20:47:35 wirawan Exp $
#
# pyqmc.sites module
#
# Wirawan Purwanto
# Created: 20110831
#

"""
pyqmc.sites

Site-specific support for pyqmc library.
"""

import re
import os
import os.path
import sys

_sites_supported = [
# Public sites
  'pmn',
# Private sites
  'wirawan0',
]

try:
  import hashlib
  def sha1digest(S):
    return hashlib.sha1(S).hexdigest()
except:
  import sha
  def sha1digest(S):
    return sha.new(S).hexdigest()

def str_grep(S, strs):
  """Returns a list of strings wherein the substring S is found."""
  return [s for s in strs if s.find(S) >= 0]

py_id_rx = re.compile(r'^[_A-Za-z][0-9_a-zA-Z]+$')
def is_python_identifier(s):
  """Verifies if this is a `valid' python identifier."""
  return py_id_rx.search(s)

def getlogin():
  try:
    return os.getlogin()
  except:
    pass
  # Fallback. This is a problem in non-login remote access, e.g.
  # using SLURM to access compute node interactively:
  try:
    return os.environ['USER']
  except:
    raise RuntimeError, "Failure getting user login info."

def _init_sites():
  """Initializes this module. Generally only needs to be called once."""
  # These variables store the contents of /etc/hosts
  # in line-by-line format (comments excluded),
  # fields-per-line list format,
  # and mines the ipv4 hostnames found.
  global _etc_hosts
  global _etc_hosts_rec
  global _etc_hosts_ipv4hosts_digest
  strip_comments = lambda S : S.split("#",1)[0].strip()
  _etc_hosts = []
  _etc_hosts_rec = []
  _etc_hosts_ipv4hosts_digest = []
  try:
    F = open("/etc/hosts", "r")
    for s in F.readlines():
      s = strip_comments(s)
      if s != "":
        _etc_hosts.append(s)
        s_fields = s.split()
        _etc_hosts_rec.append(s_fields)
        if len(s_fields) > 1:
          h1 = s_fields[0]
          if len(h1.split('.')) == 4:
            _etc_hosts_ipv4hosts_digest += [ sha1digest(hostname) for hostname in s_fields[1:] ]
    F.close()
  except:
    pass


def _autodetect_site(debug=1):
  """Performs an autodetection of the compute site where this script is
  currently running."""
  def dbg(lvl, msg):
    if debug >= lvl:
      print >> sys.stderr, msg
  for mod in _sites_supported:
    if not is_python_identifier(mod):
      # Forbid this from continuing
      dbg(1, "_autodetect_sites: Malform site name `%s'; skipped" % (mod,))
      continue
    dbg(20, "_autodetect_sites: Trying site %s" % (mod,))
    # FIXME: add try..except here:
    exec ("import pyqmc.sites.%s as _testsite" % (mod,)) in globals(), locals()
    if _testsite._detect_site():
      dbg(10, "Site detected: %s" % (_testsite.site_code))
      global site, site_config
      site = _testsite
      site_config = _testsite.site_config
      return site.site_code
  dbg(10, "_autodetect_sites: Cannot autodetect site.")
  return False


def initialize(debug=1):
  """Perform standardized initialization step.
  This needs to be called explicitly once by the main script.
  If you want to make your own hooks (e.g., custom site specification),
  you must do so before calling this routine for the first time.
  Subsequent calling should not do anything."""
  global WHEREAMI
  if "WHEREAMI" not in globals():
    where = _autodetect_site()
    if where:
      WHEREAMI = where
      return where
    else:
      return None
  return WHEREAMI


class site_config_base(object):
  """
  Generic site-specific configuration.
  Some customizable variables:
  - USER  (username on the site)
  Some root directories are also declared here:
  - GAFQMC
  - PWQMC77
  - GAMESS_CALC_ROOT
  - NWCHEM_CALC_ROOT
  - GAUSSIAN_CALC_ROOT
  - GAFQMC_CALC_ROOT
  - PWQMC77_CALC_ROOT
  - ABINIT_CALC_ROOT
  .
  """
  def __init__(self, defaults=None):
    if defaults:
      for k in defaults:
        self._defattr(k, defaults[k])
    self.init()
  def init(self):
    self._defattr('USER', getlogin())
    self.init_gafqmc()
    self.init_pwqmc()
  def _defattr(self, _attr=None, _val=None, **_kwds):
    """Sets default attribute."""
    if _attr != None and not hasattr(self, _attr):
      setattr(self, _attr, _val)
    for (k,v) in _kwds.iteritems():
      if not hasattr(self, k):
        setattr(self, k, v)
  # The following init_* routines are needed to allow
  # my legacy tools to be used.
  # They are not part of pyqmc but are often used in pyqmc-based scripts.
  # These are the most common configurations; specific sites may override
  # these defaults.
  def init_gafqmc(self):
    self._defattr(GAFQMC = os.environ['GAFQMC'])
    self._defattr( \
      GAFQMC_CALC_ROOT = self.GAFQMC,
      GAMESS_CALC_ROOT = os.path.join(self.GAFQMC, 'gamess'),
      NWCHEM_CALC_ROOT = os.path.join(self.GAFQMC, 'nwchem'),
      GAUSSIAN_CALC_ROOT = os.path.join(self.GAFQMC, 'gaussian'),
    )
  def init_pwqmc(self):
    self._defattr(PWQMC77 = os.environ['PWQMC77'])
    self._defattr( \
      PWQMC_CALC_ROOT = self.PWQMC77,
      ABINIT_CALC_ROOT = os.path.join(self.PWQMC77, 'abinit'),
    )




_init_sites()

# List of publicly 'import *'-able names:
__all__ = []
