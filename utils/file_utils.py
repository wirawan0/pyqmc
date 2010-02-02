#!/usr/bin/python
# $Id: file_utils.py,v 1.2 2010-02-02 16:11:53 wirawan Exp $
#
# pyqmc.utils.file_utils module
# File-manipulation utilities
#
# Wirawan Purwanto
# Created: 20090601
#
# Routines put here are commonly used in my own scripts.
# They are not necessarily suitable for general-purpose uses; evaluate
# your needs and see if they can them as well.
#

import bz2
import glob
import gzip
import os
try:
  import subprocess
  has_subprocess = True
except:
  has_subprocess = False


class super_file(object):
  '''"Super-file" hack wrapper for a file-like object.
  Intended to allow extra capabilities to file-like iterators such as:
  * ability to push back text lines for the subsequent next() calls.
    This is to provide some level of rewinding in parsing text files.
  * what else?
  '''
  def __init__(self, obj):
    '''Creates a super_file wrapper around the "obj" object.'''
    self.obj = obj
    self.pushback = []
  def __iter__(self):
    return self
  def close(self):
    return self.obj.close()
  def flush(self):
    return self.obj.flush()
  def next(self):
    if len(self.pushback) > 0:
      return self.pushback.pop()
    else:
      return self.obj.next()
  def push(self, s):
    self.pushback.append(s)


def open_input_file(fname, superize=0):
  if fname.endswith(".bz2"):
    fobj = bz2.BZ2File(fname, "r")
  elif fname.endswith(".gz") or fname.endswith(".Z"):
    fobj = gzip.GzipFile(fname, "r")
  elif fname.endswith(".lzma"):
    # until lzma has a "standard" python module, we use "lzma" executable:
    if has_subprocess:
      px = subprocess.Popen(("lzma", "-dc", fname), stdout=subprocess.PIPE)
      fobj = px.stdout
    else:
      px = os.popen('lzma -dc "' + fname + '"', "r")
  else:
    fobj = open(fname, "r")

  if superize:
    return super_file(fobj)
  else:
    return fobj


# Miscellaneous functions

def glob_files(filespec):
  '''Processes a glob string, or does nothing (pass-on only) if an iterable object
  (e.g. list or tuple) is already given.
  When globbing is done, the result is sorted for predictability.'''
  if getattr(filespec, "__iter__", False):
    return filespec # no re-sorting
  elif isinstance(filespec, str):
    return sorted(glob.glob(filespec))
  else:
    raise ValueError, "Don't know how to glob for an object of " + type(filespec)


def path_search(*specs, **opts):
  '''Generalized path search.
  Multiple paths can be specified for different parts of the sought filename,
  and the first file found is returned.

  Additional options:
  * pathsep="/"  -- path separator
  * filetest=os.path.isfile  -- filetest operator to be used
  * raise_error=False  -- do we want to raise an exception if the file
    is not found?
  '''
  # FIXME: this can be extremely expensive!
  xspecs = []
  xlen = []
  xstride = []
  xtot = 1
  pathsep = opts.get("pathsep", "/")
  filetest = opts.get("filetest", os.path.isfile)

  for spec in specs:
    if not getattr(spec, "__iter__", False):
      xspecs.append((spec,))
      xlen.append(1)
    else:
      xspecs.append(tuple([ x for x in spec ]))
      xlen.append(len(xspecs[-1]))
    xstride.append(xtot)
    xtot *= xlen[-1]

  for idx in xrange(xtot):
    idx0 = idx
    # Construct the filename based on the index: we reconstruct
    # the indices for all the parts given in the argument, then
    # concatenate them to get the full pathname
    s = ""
    for d in xrange(len(xspecs)-1,-1,-1):
      a = idx0 / xstride[d]
      if s == "":
        s = xspecs[d][a]
      else:
        s = xspecs[d][a] + pathsep + s
      idx0 = idx0 % xstride[d]
      #print a,
    #print s
    if filetest(s):
      return s

  if opts.get("raise_error", False):
    raise ValueError, "Cannot find file with specified combination"
  else:
    return None


def untar(archive, subdir=None, verbose=None, files=[]):
  '''Extracts a TAR archive. The destination directory can be given; otherwise
  the files are extracted to the current directory.
  Assuming GNU tar which accepts -z and -j switches.
  LZMA compression is supported via lzma program.
  '''
  opts = [ 'tar' ]
  # Python doc says: "the arguments to the child process must start with the
  # name of the command being run"

  if subdir:
    opts += [ "-C", subdir ]

  if archive.endswith(".tar.bz2") or archive.endswith(".tbz2") or archive.endswith(".tbz"):
    opts.append("-j")
  elif archive.endswith(".tar.Z") or archive.endswith(".tar.gz") or archive.endswith(".tgz"):
    opts.append("-z")
  elif archive.endswith(".tar.lzma") or archive.endswith(".tza"):
    opts.append("--use-compress-program=lzma")

  if verbose:
    for i in xrange(verbose): opts.append("-v")

  opts += [ "-xf", archive ]
  opts += files

  return os.spawnvp(os.P_WAIT, "tar", opts)


