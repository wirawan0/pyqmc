#!/usr/bin/env python
#
# pyqmc.analysis.meas_reblk.errbar1 module
#
# Created: 20150522
# Wirawan Purwanto
#
# Initial contents were imported 20150522 from Cr2_analysis_cbs.py
# (dated 20141017, CVS rev 1.143).
#

"""
pyqmc.analysis.meas_reblk.errbar1 module
Variety of tools to perform errorbar analysis of total energy.

"""

from pyqmc.analysis.cache_h5meas import H5DB_get


def reblk_errorbar(db,
                   summary=False, free_proj=True,
                   check_neg_weights=True,
                   minbeta=5.0, maxbeta=8.0, hack_reblk=False,
                   tblksizes=[2,5,10],
                   pblkcounts=[1,2,4,8,16,32,64,128,256,512,],
                   blksize=None,
                   pop_pre="none", rng=None,
                   stats_fn="reblk_errorbar.out.txt",
                   # key QMC values: must be properly given
                   nblkstep=5,
                   itv_Em=5,
                  ):
  """
  Performs combined pop-time reblocking analysis to get accurate error bar estimate.

  This was originally subroutine `reblk_Cr2_errbar` in Cr2 project,
  developed to aid errorbar analysis of a free projection run.

  Default parameters were from Cr2 project (~spring 2014).
  In those FP runs, the measurements were done every beta=0.1
  (and the nblkstep also corresponds to beta=0.1).
  You better adjust these to your own case!
  """
  from numpy import average, product, std, sqrt, sum
  import pyqmc.analysis.meas_reblk.H5meas_reblk as h5m

  h5db = H5DB_get(db)

  # opts
  pass
  # key values
  nblkstep = 5
  itv_Em = 5
  nblk_Em = nblkstep // itv_Em

  def matchme(m):
    return (minbeta is None or m['beta'] >= minbeta) \
       and (maxbeta is None or m['beta'] <= maxbeta)

  #free_proj = True

  if not tblksizes:
    tblksizes = [nblk_Em]

  if not pblkcounts:
    pblkcounts = [1]

  h5blk = h5m.h5meas_reblk(h5db, free_proj=free_proj,
               tblksizes=tblksizes,
               pblkcounts=pblkcounts,
               # Filter for records:
               match=matchme, stage=2,
               blksize=blksize, pop_pre=pop_pre, rng=rng,
               )

  if check_neg_weights:
    # TODO
    pass
    #h5meas_reblk_check_negative_weight(h5blk)

  if summary:
    if False:
      print "Raw data:"
      for r in h5blk1.raw_meas:
        beta, E, w = r[1], r[2], r[3]
        print beta, w, E

    if isinstance(summary, tuple):
      pblk, tblk = summary
    else:
      pblk = pblkcounts[0]
      tblk = tblksizes[0]
    print "Showing summary for pblk,tblk = %s,%s" % (pblk, tblk)
    rawblk = h5blk.raw_blocks[pblk, tblk]
    #rawblk_ndata = rawblk.shape[1] # num of data across time axis
    rawblk_ndata = product(rawblk.shape) # num of data across ALL axes
    # raw.out-style reblocked data:
    print "Reblocked data:"
    for (r,Xw,w) in zip(h5blk.raw_meas[tblk-1::tblk], rawblk['Xw'].T, rawblk['w'].T):
      print "%3.1f %.12f %.13f" % (r['beta'], sum(w), sum(Xw)/sum(w))

    #for (r,Xw,w) in zip(h5blk1.raw_meas[tblk-1::tblk], rawblk['Xw'][0], rawblk['w'][0]):
    #  print "%3.1f %.12f %.13f" % (r['beta'], w, Xw/w)

    print "Grand average:", sum(rawblk['Xw']) / sum(rawblk['w'])
    print "Block average:", average(rawblk['Xw'] / rawblk['w'])
    print "Block error:  ", std(rawblk['Xw'] / rawblk['w'], ddof=1) / sqrt(rawblk_ndata)
    ANS = """
    """
  if stats_fn:
    h5m.h5meas_dump_reblk_stats(h5blk, stats_fn)

  return h5blk

