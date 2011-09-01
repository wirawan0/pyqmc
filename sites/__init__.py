#!/usr/bin/python
# $Id: __init__.py,v 1.2 2011-09-01 20:10:11 wirawan Exp $
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

_sites_supported = [
  'pmn',
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

_init_sites()

# List of publicly 'import *'-able names:
__all__ = []
