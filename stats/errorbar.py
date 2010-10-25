# $Id: errorbar.py,v 1.3 2010-10-25 14:46:58 wirawan Exp $
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
    # Standard errorbar matcher (errorbar in parantheses), in the form of a tuple
    # The first element of the matcher is a regexp matcher object, and the
    # second element is a function to extract the (value,error) tuple from the
    # successful matching process.
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
    # Later additional matchers can be added here
    R.init = True

  @staticmethod
  def aux():
    R = regexp__aux
    if not R.init: R.initialize()
    return R

  @staticmethod
  def match(matcher, Str, rslt, flatten=False):
    '''Matches the string `Str' against the errorbar regexp pattern in `matcher[0]'.
    If it matches, the value and error are added to the `rslt' list.
    Depending on whether `flatten' is True or not, it is added as a tuple or as
    two elements, respectively, into the `rslt' list.'''
    # Note: matcher is an object like R.errbar above.
    m = matcher[0].match(Str)
    if m:
      if flatten:
        rslt.extend(matcher[1](m))
      else:
        rslt.append(matcher[1](m))
      return True
    else:
      return False

def expand(obj, convert_float=False, flatten=False):
  '''Expands compressed errorbar notation to two consecutive numbers
  (returned as a tuple).

  Input: The input can be a string, or a list of strings.

  Output:
  The list element that has the format of "VAL(ERR)" (or its scientific
  format twist) will be expanded into two numbers in the output list.
  All other elements will be passed "as is" to the output.
  Optionally, the non-float items can be force-converted to floats if
  the convert_float is set to True.
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
      if (rgx.match(rgx.errbar, o, rslt, flatten)):
        #print "match: errbar"
        pass
      elif convert_float:
        # Convert to float right away, store into the `rslt' list
        rslt.append(float(o))
        #rslt.append( (float(o),) )
      else:
        # Unless otherwise requested, the object will not be converted
        # to float:
        rslt.append(o)
  return rslt

class float_decomposer(object):
  """Floating-point decomposition.
  We are assuming IEEE double precision here."""

  def __init__(self, val):
    self.val = val
    V = "%+.16g" % val
    self.sign = V[0]
    self.digits = V[1] + V[3:3+16]
    self.exp = int(V[18:])

  set = __init__


class errorbar_compressor(object):
  """Compressor for errorbar string."""
  def __init__(self):
    self.errdigits = 2
  def __call__(self, val, err, **args):
    errdigits = args.get("errdigits", self.errdigits)


