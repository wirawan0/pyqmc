# $Id: result_base.py,v 1.1 2010-10-25 15:29:42 wirawan Exp $
#
# pyqmc.results.result_base.py
# Basic tool for storing the result of a calculation.
#
# Wirawan Purwanto
# Created: 20101025
#
# The class result_base is very similar to
#   wpylib.params.params_flat.Parameters
# except that it does not contain multiple dict-like objects, thus is
# much simpler.
#

class result_base(dict):
  '''Structure to represent metadata or structured result.
  No particular structure is assumed.
  Results are fetchable by either X.member or X['member'] syntax, and these
  results are typically meant to be read-only.

  CAVEATS
  * Note: dict method method names are left intact.
    Please be aware when programming this thing!
  * __setattr__ is not set.
    Result-related attributes are supposed to be read-only.
  * As a consequence, adding metadata or result *must* be done in the
    dict way, e.g.
        X['nblocks'] = 32
    instead of
        X.nblocks = 32.
  '''
  def __init__(self, src=None):
    if isinstance(src, dict):
      self.clear()
      self.update(src)
    if isinstance(src, str):
      # WARNING: Awaiting future definition of parse_text_file_().
      # This must be specified in the derived class.
      self.parse_text_file_(src)
      self.filename_ = src
    else:
      pass
  def __getattr__(self, key):
    try:
      return self[key]
    except:
      return dict.__getattribute__(self, key)
  def __str__(self):
    return "<" +self.__module__ + "." + self.__class__.__name__ + \
      " object (" + dict.__str__(self) + ")>"


# Test programs
if __name__ == "__main__":

  def _test_result_base1():
    x = result_base()
    print x.keys()
    print x.values()
    x['nblocks'] = 32
    print x

    y = result_base({'lack': 352, 'Etrial': -32.764})
    print y.keys()
    print y.values()
    y['nblocks'] = 32
    print y
    print y.Etrial


  _test_result_base1()
