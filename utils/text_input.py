#!/usr/bin/python
# $Id: text_input.py,v 1.1 2009-06-12 15:06:50 wirawan Exp $
#
# pyqmc.utils.text_input module
# Quick-n-dirty text input utilities
#
# Wirawan Purwanto
# Created: 20090601
#
# Routines put here are commonly used in my own scripts.
# They are not necessarily suitable for general-purpose uses; evaluate
# your needs and see if they can them as well.
#

import re

import numpy

from pyqmc.utils.file_utils import open_input_file


class text_input(object):
  '''Text input reader with support for UNIX-style comment marker (#) and
  standard field separation (tabs and whitespaces).
  Used for quick and dirty data reading (iterating only once in forward
  direction without the need of rewinding or skipping).
  This object can be treated like an input file, e.g. used as an iterator,
  etc.

  To support more fancy options (e.g., rewinding), use "superize=1" when
  creating the instance.'''

  def __init__(self, fname, **opts):
    if opts.get("superize", 0):
      open_opts = { "superize" : opts["superize"] }
      del opts["superize"]
    else:
      open_opts = {}
    self.file = open_input_file(fname, **open_opts)
    # field_filtering_proc field can be used to filter unwanted fields, or do
    # some additional transformations before final feed to the main iteration.
    self.field_filtering_proc = lambda flds : flds
    if len(opts) > 0:
      self.set_options(**opts)

  def __del__(self):
    self.file.close()

  def __iter__(self):
    return self

  """
  def next(self):
    while True:
      L = self.file.next()
      F = self.field_filtering_proc(L.split("#")[0].split())
      if len(F) > 0:
        return F
  """

  def next_rec(self):
    '''Yields the next record, which is already separated into fields.'''
    while True:
      L = self.file.next()
      F = self.field_filtering_proc(L.split("#")[0].split())
      if len(F) > 0:
        return F

  def next_line(self):
    '''Yields the next line, which is already separated into fields.'''
    while True:
      L = self.file.next()
      F = self.field_filtering_proc(L.split("#")[0].rstrip())
      if len(F) > 0:
        return F

  # Do NOT touch the "next" field below unless you know what you're doing:
  next = next_line

  def seek_text(self, regex=None, match=None):
    '''Seeks the file until a particular piece text is encountered.
    We ignore all comments.'''

    if regex:
      if isinstance(regex, str):
        Regexp = re.compile(regex)
      else:
        Regexp = regex
      match_proc = lambda x: Regexp.search(x)
    else:
      match_proc = match

    while True:
      L = self.next_line()
      if match_proc(L):
        return L


  def read_floats(self, *cols, **kwd):
    """Quickly reads a set of floats from a text file.
    Returns a numpy array of the values in double precision.

    Example usage:
      >>> arr = text_input("/tmp/file.txt").read_floats(0, 2, 3)
    to read columns 1, 3, and 4 of the text file /tmp/file.txt, while disregarding
    comments.
    """
    # float_fields extracts the desired columns and converts them to floats
    float_fields = lambda vals : [ float(vals[col]) for col in cols ]
    if "maxcount" in kwd:
      rslt = [ float_fields(vals.split()) for (c,vals) in zip(xrange(kwd['maxcount']),self) ]
    else:
      rslt = [ float_fields(vals.split()) for vals in self ]
    # finally convert them to a numpy ndarray:
    return numpy.array(rslt)

  def read_items(self, *col_desc, **kwd):
    """Quickly reads a set of items from records of whitespace-separated fields
    in a text file.
    Returns a structured numpy array of the values read.

    Example usage:

      >>> arr = text_input("/tmp/file.txt").read_items(0, (2, int), (3, "S10", "Atom"))

    reads columns 1 (as floats, by default), 3 (as integers), and 4 (as strings of
    max length of 10, which field is named "Atom") from the text file /tmp/file.txt,
    while disregarding comments.

    If the tuple contains the third field, it is used as the name of the field;
    otherwise the fields are named f0, f1, f2, ....


    TODO: Needs ability to read in complex data.
    """
    deftype = kwd.get("deftype", float)

    # float_fields extracts the desired columns and converts them to floats
    flds = []
    cols = []
    for (i,c) in zip(xrange(len(col_desc)), col_desc):
      if type(c) == int:
        cols.append(c)
        flds.append(('f' + str(i), deftype))
      elif len(c) == 1:
        cols.append(c[0])
        flds.append(('f' + str(i), deftype))
      elif len(c) == 2:
        cols.append(c[0])
        flds.append(('f' + str(i), c[1]))
      elif len(c) == 3:
        cols.append(c[0])
        flds.append((c[2], c[1]))

    #print cols
    #print flds
    get_fields = lambda vals : tuple([ vals[col] for col in cols ])
    if "maxcount" in kwd:
      print "hello"
      rslt = [ get_fields(vals.split()) for (c,vals) in zip(xrange(kwd['maxcount']),self) ]
    else:
      rslt = [ get_fields(vals.split()) for vals in self ]
    #print rslt
    # finally convert them to a numpy ndarray:
    return numpy.array(rslt, dtype=flds)

  # Sets fancy options
  def set_options(self, **opts):
    for (o,v) in opts.iteritems():
      if o == "expand_errorbar":
        self.expand_errorbar(v)
      else:
        raise "ValueError", "Invalid option: %s" % (o,)
    return self

  # Option for errorbar expansion:
  def expand_errorbar(self, v=True):
    '''Enables or disables errorbar expansion.'''
    if v:
      self.opt_expand_errorbar = True
      self.field_filtering_proc = self.expand_errorbar_hook
    else:
      self.opt_expand_errorbar = False
      self.field_filtering_proc = lambda flds : flds
    return self

  def expand_errorbar_hook(self, F):
    # A hook for field_filtering_proc for expanding errorbars:
    from pyqmc.stats.errorbar import expand
    return expand(F, flatten=True)
