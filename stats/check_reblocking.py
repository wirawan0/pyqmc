#!/usr/bin/python
# $Id: check_reblocking.py,v 1.2 2009-02-05 05:56:36 wirawan Exp $
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
import numpy

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
    ('ndata', '=i4'),
    ('blksize', '=i4'),
    ('reblk_ndata', '=i4'),
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

  Option plot_out can be true for default (dumb text) plotting, or an existing
  Gnuplot.Gnuplot() instance.
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
    E_blk_stddev = reblk_E.std() * sqrt(reblk_nblk / (reblk_nblk-1.0))
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
