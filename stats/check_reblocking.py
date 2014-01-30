#!/usr/bin/python
# $Id: check_reblocking.py,v 1.3 2010-02-24 15:20:53 wirawan Exp $
#
# check_reblocking.py
# Python implementation of my famous "check_reblocking" tool.
#
# Wirawan Purwanto
# Created: 20081022
#
# The use of numpy is all that takes to speed up the calculation
# tremendously! Implementing as a pure perl (or python) script is
# very portable (having virtually no dependency) but it runs very
# slow for large datasets.
#
import sys
import numpy

from wpylib.db.result_base import result_base

NDATA = 0
BLKSIZE = 1
REBLK_NBLK = 2
UPTOT = 3
DOWNTOT = 4
E_BLK_AVG = 5
E_BLK_STDDEV = 6
E_BLK_ERR = 7

check_reblocking_type = \
    numpy.dtype([
    ('ndata', '=i8'),
    ('blksize', '=i8'),
    ('reblk_ndata', '=i8'),
    ('uptot', '=f8'),
    ('downtot', '=f8'),
    ('Eavg', '=f8'),
    ('Estddev', '=f8'),
    ('Eerr', '=f8'),
    ])

def check_reblocking(wtwlkr, El, reblk_sizes, print_out = False, \
  print_timing = False, plot_out = False, format = 1): # , El_caps = None):
  '''Performs reblocking check.

  El_caps is an optional 2-number tuple which specifies the minimum and maximum
  energies allowed for the measurement.

  Returns a list of 8-component tuples, each of which containing

  (ndata, blksize, reblk_nblk, uptot, downtot, E_blk_avg, E_blk_stddev, E_blk_err)

  The list entries correspond to the block sizes in reblk_sizes input array.

  If print_out==True, it will print the statistics in text tabular form,
  containing the following fields:

      ndata, blksize, reblk_nblk, E_blk_avg, E_blk_stddev, E_blk_err

  Option plot_out can be a true value, which will display the plot in a
  text format (dumb terminal), or an existing Gnuplot.Gnuplot() instance.
  '''
  from numpy import array, sqrt
  from time import clock
  import Gnuplot

  if print_timing: clk1 = clock()
  # The input data must be 1-D arrays (or list)!
  wtwlkr = array(wtwlkr) # data["wtwlkr"])
  #print type(wtwlkr)
  El_w = array(El) * wtwlkr
  ndata = El_w.shape[0]

  '''# --- This capping should not have been done here!
  if El_caps != None:
    (El_min, El_max) = El_caps
    ncap_lo = 0
    ncap_hi = 0
    print "applying E_l caps from %g to %g" % (El_min, El_max)
    for (i, E) in zip(xrange(ndata), El):
      if E < El_min:
        ncap_lo += 1
        El_w[i] = El_min * wtwlkr[i]
      elif E > El_max:
        ncap_hi += 1
        El_w[i] = El_max * wtwlkr[i]
    if ncap_hi > 0 or ncap_lo > 0:
      print "Total %d elements capped: %d too high, %d too low" % \
        (ncap_hi + ncap_lo, ncap_hi, ncap_lo)
  '''

  Rslt = []
  errbars = []

  for blksize in reblk_sizes:
    # rblk_nblk = number of new blocks as a result of reblocking
    # rblk_ndata = number of data points usable as a result of reblocking
    reblk_nblk = int(ndata / blksize)
    reblk_ndata = reblk_nblk * blksize
    # Reshape the arrays for quick averaging:
    reblk_wt = wtwlkr[0:reblk_ndata].reshape(reblk_nblk, blksize).sum(1)
    reblk_El_w = El_w[0:reblk_ndata].reshape(reblk_nblk, blksize).sum(1)
    # This is the reblocked energies: array[0:reblk_nblk]
    reblk_E = reblk_El_w / reblk_wt
    #print reblk_E[0:20]

    E_blk_avg = reblk_E.mean()
    # We must use unbiased variance & stddev:
    if reblk_nblk > 1:
      E_blk_stddev = reblk_E.std() * sqrt(reblk_nblk / (reblk_nblk-1.0))
    else:
      E_blk_stddev = 0
    E_blk_var = E_blk_stddev**2
    # reblk_E.std()**2 * reblk_nblk / (reblk_nblk-1)
    E_blk_err = E_blk_stddev / sqrt(reblk_nblk)
    # "Grand" averaging:
    uptot = sum(reblk_El_w)
    downtot = sum(reblk_wt)
    E_grand_avg = uptot / downtot

    # Return the reduction results in a tuple
    rslt = (ndata, blksize, reblk_nblk, uptot, downtot, E_blk_avg, \
            E_blk_stddev, E_blk_err)
    #rslt = (ndata, blksize, reblk_nblk, reblk_ndata, E_blk_avg, E_grand_avg, \
    #        E_blk_stddev, E_blk_err)
    Rslt.append(rslt)
    errbars.append(E_blk_err)

  if plot_out:
    ddd = dir(plot_out)
    if 'splot' in ddd and 'gnuplot' in ddd:  # a Gnuplot.Gnuplot() instance?
      g = plot_out
    else:
      g = Gnuplot.Gnuplot()
      g("set term dumb")
    d = Gnuplot.Data(reblk_sizes, errbars, with_ = "linespoints")
    g.plot(d)

  if print_out:
    for rslt in Rslt:
      print "%7d %7d %7d %17.10g %17.10g %17.10g" % \
        (rslt[0], rslt[1], rslt[2], rslt[5], rslt[6]**2, rslt[7])

  if print_timing:
    clk2 = clock()
    print "total time = %.4f secs" % (clk2 - clk1)

  if format == 0:  # legacy format: no array-izing
    return Rslt
  else:
    return numpy.array(Rslt, dtype=check_reblocking_type)



class weighted_samples(object):
  """Representation of set of measurement samples with associated
  weights.
  This OO interface is an alternative to, and richer implementation of,
  the check_reblocking subroutine.
  """
  def init(self, X, w, copy=True):
    self.X = numpy.array(X, copy=copy)
    self.w = numpy.array(w, copy=copy)
    self.X_w = self.X * self.w
    assert self.X.shape == self.w.shape
    assert len(self.X.shape) == 1  # enforce 1D data for now.
    self.ndata = self.X.shape[0]

  def reblock(self, blksize, shift=0, transpose=False, save=False, stats="std"):
    """Reblock the samples according to block size `blksize'.

    The `shift' argument, if given, must be one of the following:
    - a real number, 0 <= shift < 1, representing the precentage of
      remaining (unblocked) part of the sample array.
    - an integer >= 1, representing the number of data points skipped
      at the beginning of the array.

    The `save' argument, if evaluates to Boolean True, saves the calculated
    values & intermediate arrays to _reblk_rec array.

    The `stat' argument specifies the extra statistical processings done
    on the reblocked data:
    - 'std' computes variance, standard deviation, and standard error.

    """
    from numpy import average, std, sqrt
    ndata = self.ndata
    assert blksize > 0
    assert blksize <= ndata
    assert shift >= 0
    # rblk_nblk = number of blocks as a result of reblocking
    # rblk_ndata = number of effective (usable) data points as
    # a result of reblocking
    reblk_nblk = int(ndata // blksize)
    reblk_ndata = reblk_nblk * blksize
    shift_max = ndata - reblk_ndata

    if 0 < shift < 1.0:
      shift = int(round(shift * shift_max))
    else:
      if not isinstance(shift, (int,long)):
        shift = int(round(shift))
    assert shift <= shift_max

    # Reshape the arrays for quick averaging:
    if not transpose:
      blk_w = self.w[shift : reblk_ndata+shift].reshape(reblk_nblk, blksize)
      blk_Xw = self.X_w[shift : reblk_ndata+shift].reshape(reblk_nblk, blksize)
    else:
      blk_w = self.w[shift : reblk_ndata+shift].reshape(blksize, reblk_nblk).T
      blk_Xw = self.X_w[shift : reblk_ndata+shift].reshape(blksize, reblk_nblk).T

    reblk_w = blk_w.sum(1)
    reblk_X_w = blk_Xw.sum(1)
    # This is the reblocked averages: array[0:reblk_nblk]
    reblk_X = reblk_X_w / reblk_w

    if save:
      #L = locals()
      rec = result_base(
        raw_ndata=ndata,
        blksize=blksize,
        shift=shift,
        shift_max=shift_max,
        reblk_nblk=reblk_nblk,
        reblk_ndata=reblk_ndata,
        blk_w=blk_w,
        blk_Xw=blk_Xw,
        reblk_w=reblk_w,
        reblk_X_w=reblk_X_w,
        reblk_X=reblk_X,
      )
      self._reblk_rec = rec
      if stats in ('std',):
        rec['mean'] = average(rec.reblk_X)
        rec['std'] = std(rec.reblk_X, ddof=1)
        rec['var'] = rec.std**2
        rec['err'] = rec.std / sqrt(rec.reblk_nblk)

    return (reblk_X, reblk_w)

  def shuffle(self, rng=None):
    """Shuffles the data arrays in-place, in a coherent manner (all fields:
    X, w, X_w).
    """
    if rng == None:
      rng = numpy.random
    state = rng.get_state()
    rng.shuffle(self.X)
    rng.set_state(state)
    rng.shuffle(self.w)
    rng.set_state(state)
    rng.shuffle(self.X_w)



class reblocking_result(object):
  """An array to hold results from reblocking tests.
  Useful for collecting results from reblocking experiment in a
  numpy-style tabular form.

  Limitation: the number of records must be known ahead of time.
  """
  dtype_map = {
    'raw_ndata': int,
    'blksize': int,
    'reblk_nblk': int,
    'shift': int,
    'shuffle': int,  # shuffling index; -1 means no shuffling
    'mean': float,
    'var': float,
    'std': float,
    'err': float,
  }
  def __init__(self, numrecs, Dtype):
    """Initializes a new reblocking result table.
    Input:

    - numrecs (int): number of records to add.
    - Dtype: an array of two-tuples defining the fields of the table as well
      as the output formatting of the row in a tabular text form.

    Example:

    >>> Rec = reblocking_result(number_of_rows,
          [('raw_ndata',    '%9d'),
           ('blksize',    '  %8d'),
           ('reblk_nblk', ' %10d'),
           ('shift',      '  %8d'),
           ('mean',    '   %16.10f'),
           ('var',      '  %#16.8f'),
           ('err',      '  %12.10f')]
        )

    The datatype and field names are hardcoded--see dtype_map attribute
    of this class.
    """
    dm = self.dtype_map
    self.fields = [ fld for (fld, fmt) in Dtype ]
    self.fmt = "".join([ fmt for (fld, fmt) in Dtype ]) + "\n"
    self.dtype = [ (f, dm[f]) for f in self.fields ]
    self.nrows = 0
    self.arr = numpy.zeros(numrecs, dtype=self.dtype)
    #return self.arr

  def append(self, rec):
    """
    Input: rec is the data structure returned by weighted_samples.reblock()
    method, as `weighted_samples._reblk_rec' field.
    """
    dest = self.arr[self.nrows:self.nrows+1]
    for f in self.fields:
      dest[f][0] = rec[f]
    self.nrows += 1

  def print1(self, rec, out=sys.stdout):
    """Prints out a single record.

    Input: rec is the data structure returned by weighted_samples.reblock()
    method, as `weighted_samples._reblk_rec' field.
    """
    rec1 = tuple( rec[f] for f in self.fields )
    out.write(self.fmt % rec1)

  def table_header(self):
    """Prints out the table header for text tabular printout."""
    from wpylib.text_tools import str_fmt_heading
    fmt = str_fmt_heading(self.fmt)
    return fmt % tuple(self.fields)

  def __getitem__(self, i):
    """Returns the i-th record of the table."""
    return self.arr[i]



def make_sorted_unique_int(series):
  """Given a list of numbers, round them to nearest integers, sort them,
  and return only the unique values.
  """
  s_rounded = numpy.array(numpy.round(series),dtype=int)
  sorted_list = list(set(s_rounded))
  sorted_list.sort()
  return numpy.array(sorted_list)


def reblock_sizes(ndata, nblksizes=50, minrblks=20):
  """Generate unique reblock block sizes on a logarithmic scale (actually,
  an exponential scale) to study reblocking sizes.

  Input:
  - ndata = number of raw data points
  - minrblks = minimum number of reblocked data points desired at the largest
    block size.
    The largest block size is thus approximately (ndata / minrblks)

  Example:
  >>> reblock_sizes(60000,nblksizes=50,minrblks=20)
  array([   1,    2,    3,    4,    5,    6,    7,    8,   10,   12,   14,   16,   19,   22,   26,   31,   36,   43,   50,   59,   70,   82,   97,  114,  135,  158,  187,  220,  259,  305,  359,  422,  497,  585,  689,  812,  956, 1126, 1325, 1561, 1838, 2164, 2548, 3000])

  """
  maxblksize = float(ndata) / minrblks
  blks = numpy.exp(numpy.linspace(0, numpy.log(maxblksize), nblksizes))
  return make_sorted_unique_int(blks)

