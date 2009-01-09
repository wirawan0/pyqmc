# $Id: errorbar.py,v 1.1 2009-01-09 22:06:16 wirawan Exp $
#
# Module pyqmc.stats.errorbar
# Errorbar text handling for Python
#
# Created: 20081118
# Wirawan Purwanto
#

import re

class regexp__aux(object):
  '''Auxiliary objects for routines below. Compiled into regexp objects
  for speed.'''
  # CHECKME FIXME: This class is NOT originally designed to be multithread-safe
  # since regex objects are shared between all threads.
  # To be thread-safe, the regex objects or aux classes may need to be
  # instantiated separately for each thread instance, or whatever else way
  # possible.
  init = False

  @staticmethod
  def initialize():
    R = regexp__aux
    # Each of the regex stuff below is a 2-tuple containing the regex
    # object and its associated function to extract the value/errbar pair
    # (as a tuple) from the re.Match object.
    #
    # Standard errorbar matcher (errorbar in parantheses)
    R.errbar = \
      (
        re.compile(r"([-+]?\d+)"    # leading digits with optional sign
                   r"((?:\.\d*)?)"  # optional digits after decimal point
                   r"\((\d+\.?\d*[Ee]?[-+]?\d*)\)" # errorbar in paranthesis
                   r"((?:[Ee][-+]?\d+)?)$"),
        lambda M: ( float(M.group(1) + M.group(2) + M.group(4)), # value = simple composition
                    # for errorbar:
                    # - must multiply with value's exponent
                    # - and account for decimal digits (if any)
                    float(M.group(3)) * float("1"+M.group(4)) * 10**(-max(len(M.group(2))-1, 0)) )
      )
    # Later additional matchers can be added

    R.init = True

  @staticmethod
  def aux():
    R = regexp__aux
    if not R.init: R.initialize()
    return R

  @staticmethod
  def match(matcher, Str, rslt):
    m = matcher[0].match(Str)
    if m:
      rslt.append(matcher[1](m))
      return True
    else:
      return False

def expand(obj):
  '''Expands compressed errorbar notation to two consecutive numbers
  (returned as a tuple).

  Input: The input can be a string, or a list of strings.

  Output:
  The list element that has the format of "VAL(ERR)" (or its scientific
  format twist) will be expanded into two numbers in the output list.
  All other elements will be converted to floats and passed "as is" to the
  output.
  '''
  if getattr(obj, "__iter__", False):
    iterable_inp = True
    objs = obj
  else:
    iterable_inp = False
    objs = ( obj, )

  rgx = regexp__aux.aux()

  rslt = []
  for o in objs:
    t = type(o)
    if t == int  or  t == float  or  t == complex:
      rslt.append(o)
    else:
      # Assume a string!
      o = o.strip()
      #m = rgx.errbar.match(o)
      if (rgx.match(rgx.errbar, o, rslt)):
        #print "match: errbar"
        pass
      else:
        # Convert to float right away; store as 1-tuple instead of
        # 2-tuple
        rslt.append( (float(o),) )
  return rslt

