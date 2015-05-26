#!/usr/bin/env python
#
# pyqmc.analysis.meas_reblk.H5meas_reblk module
#
# Created: 20150522
# Wirawan Purwanto
#
# Initial contents were imported 20150522 from Cr2_analysis_cbs.py
# (dated 20141017, CVS rev 1.143).
#


"""
pyqmc.analysis.meas_reblk._h5_meas_reblk
QMC energy measurement reblocking tools

This is a more powerful and in-depth approach to post-runtime
QMC errorbar analysis than the raw check_reblocking approach.


Note:
The dataset must have been converted to HDF5 archival format before
processed.


"""

import numpy

from wpylib.params import struct

class meas_reblk_result(struct):
  """A class to hold (total energy measurement) reblocking results.
  """

# Imported 20150522 from Cr2_analysis_cbs.py
# (dated 20141017, CVS rev 1.143).

def h5meas_reblk(meas_db, free_proj=False,
                 tblksizes=None,
                 pblkcounts=None,
                 # Filter for records:
                 match=None, # (lambda m : True),
                 #minbeta=0, maxbeta=None,  ## redundant! use match routine!
                 # advanced loader tweaks, memory optimization
                 blksize=None, stage=None,
                 pop_pre="none", rng=None,
                 #adv_reblk=True, reblk_sizes=None,
                 reblk_save_logs=True,
                 reblk_log_root_dir=".",
                 debug=100,
                ):
  """Analyze errorbar from hdf5 measurement record (for time-population
  combined reblocking).
  The primary input, `meas_db`, is a python object that represents
  the raw QMC energy measurements in HDF5 format, which is of one of
  these classes:

      pyqmc.results.gafqmc_meas.meas_hdf5
      pyqmc.results.pwqmc_meas.meas_hdf5

  If the data is from free projection run,m then set `free_proj` to True.
  Otherwise, it is assumed to be a phaseless-constrained run.

  In this analysis, the records (which can be thought roughly as
  a two-dimensional array of energy values, where the two axes are
  the imaginary time and walker (a.k.a. population) axes.

  The dataset is reblocked in imaginary-time axis in the usual manner:
  by giving the sequence of block sizes in `tblksizes` input (1-D) array.
  The dataset is reblocked in walker (population) axis in a different way:
  we divide up the measurement data into several roughly equal partitions;
  the partition sizes to use is given in the `pblkcounts` input (1-D) array.
  Because the number of walkers fluctuates during the simulation,
  the length in the population axis is not fixed, thus the equal partitioning
  is an approximation.

  NOTE: The pop-axis reblocking is slightly different from reblock_meas()
  in Kpt-MnO-AFM2.py because instead of "throwing out" excess data which
  do not make it when we "rectangularly" reshape the array, we include
  all of them by allowing some of the partition to be slightly larger
  (excess by 1 walker from the other).
  Statistically in the limit of very large number of walkers, the two should
  be equal.

  The loading of the measurement data into memory can be blocked
  to conserve memory usage during the processing.
  Use `blksize` parameter for that purpose; then only `blksize` number of rows
  are loaded at a time.
  If `blksize` is a literal string 'block', then the data will be loaded
  block-by-block, as defined in QMC measurement block.

  The data can be filtered only for a given `stage` (which is what
  you want if blksize=='block'):
    0 == equilibration,
    1 == growth,
    2 == measurement phase

  There are several built-in measurement data preprocessing methods
  that can be applied right after the data is read from the file but
  before the block-averaging is done. This is determined by the string
  value of `pop_pre` argument, which can be one of these:

    'none' or 'contiguous'
    'shuffle'
    'coherent-shuffle'
    'buggy-contiguous'  (not recommended; can be problematic)
    'buggy-transpose'   (not recommended; can be problematic)

  For shuffle operations, the `rng` argument can be given (and is like
  the numpy-style `numpy.random.RandomState` RNG object) so the shuffle
  order is reproducible. This is the recommended approach.


  DINOSAUR TRACKS

  Please do not use arguments `reblk_save_logs` and `reblk_log_root_dir`.

  """
  from sys import stderr
  from itertools import product
  from wpylib.math.stats import jackknife1
  from wpylib.math.stats.jackknife1 import jk_generate_averages, jk_stats_aa
  from wpylib.text_tools import str_fmt_heading
  from pyqmc.stats.check_reblocking import \
    weighted_samples, reblocking_result, \
    reblock_sizes
  #from reblocking import Reblock_contiguous, Reblock_interleaving, Reblock_shuffled
  global ws

  def dbg(lvl, msg):
    if debug >= lvl:
      stderr.write(msg)

  match_params = {} # TBD
  pblkcounts = numpy.asarray(pblkcounts, dtype=int)
  tblksizes = numpy.asarray(tblksizes, dtype=int)
  # callback funcs
  def match_rows(m):
    row = m['row']
    return row >= match_params['start'] and row < match_params['stop'] \
      and (match is None or match(m))
      #and (minbeta is None or m['beta'] >= minbeta) \
      #and (maxbeta is None or m['beta'] <= maxbeta)
  def match_block(m):
    blk = m['iblk']
    return blk == match_params['start'] \
      and (match is None or match(m))

  raw_dset = meas_db.raw

  dbg(1, "h5meas_reblk: db=%s\n" % (meas_db.filename()))
  dbg(2, "h5meas_reblk: pop_pre=%s\n" % (pop_pre))

  # Figure out how many datasets that match the "match" filter
  # requirements in the raw dataset
  nrows = 0
  match_row_begin = None
  match_row_end = None
  meta = raw_dset.meta[:]
  if match is None:
    match1 = lambda m : True
  else:
    match1 = match

  num_wlkrs = []
  if stage is None:
    dbg(100, "matching: route 3\n")
    for (irow,m) in enumerate(meta):
      if m['count'] > 0 and match1(m):
        if match_row_begin is None: match_row_begin = irow
        match_row_end = irow
        nrows += 1
        num_wlkrs.append(m['count'])
  else:
    dbg(100, "matching: route 4\n")
    for (irow,m) in enumerate(meta):
      if m['count'] > 0 and m['stage'] == stage and match1(m):
        if match_row_begin is None: match_row_begin = irow
        match_row_end = irow
        nrows += 1
        num_wlkrs.append(m['count'])
  dbg(5, "Matching rows: %d (%d .. %d)\n" % (nrows, match_row_begin, match_row_end))
  num_wlkrs = numpy.array(num_wlkrs)
  dbg(5, "Num walkers: [%d .. %d]  avg=%.2f stddev=%.2f\n" \
         % (num_wlkrs.min(), num_wlkrs.max(), num_wlkrs.mean(), num_wlkrs.std(ddof=1)))
  num_wlkrs_max = num_wlkrs.max()

  if blksize is None:
    # this can be >= len of actual measurement data, but that's ok.
    blksize = nrows # raw_dset.row_count()

  if blksize == 'block':
    match_params['start'] = 0
    match_params['skip'] = 1
    match_proc = match_block
  else:
    match_params['start'] = match_row_begin
    match_params['stop'] = match_row_begin + blksize
    match_params['skip'] = blksize
    match_proc = match_rows

  if reblk_save_logs:
    def logfile(*fn):
      return path_prep(reblk_log_root_dir, *fn)
  else:
    def logfile(*fn):
      return None

  raw_blk_dtype = [ ('Xw', float), ('w', float) ]
  raw_planks = dict( (pb, numpy.zeros((pb, nrows), dtype=raw_blk_dtype)) for pb in pblkcounts )
  raw_blocks = dict( ((pb,tb), numpy.zeros((pb, nrows // tb), dtype=raw_blk_dtype)) for (pb,tb) in product(pblkcounts,tblksizes) )
  #raw_blocks = ()

  # Some basic statistics on raw measurement (per time snapshot)
  raw_meas_dtype = [
    ('count', int), ('beta', float),
    ('X', float), ('w', float),
    ('Xjk', float), ('Xjk_err', float),
  ]
  raw_meas = numpy.zeros((nrows,), dtype=raw_meas_dtype)

  #Fsum = text_output(logfile("summary-reblocking.dat.txt"), flush=True)

  if rng is None:
    rng = numpy.random

  i = 0
  ws = weighted_samples()
  has_header = False
  rslt = meas_reblk_result()
  rslt.raw_numrows = nrows
  rslt.pblkcounts = pblkcounts
  rslt.tblksizes = tblksizes
  rslt.raw_meas = raw_meas
  rslt.raw_planks = raw_planks
  rslt.raw_blocks = raw_blocks
  rslt.db_filename = meas_db.filename()

  act_row_begin = None
  act_row_end = None

  __first = False   ## ONLY ENABLE FOR SEVERE DEBUGGING

  if pop_pre == 'coherent-shuffle':
    """Determines the overall order that essentially is independent on the
    number of walkers, but approximately preserves the order of the walkers
    from one time snapshot to another."""
    #wlk_order = numpy.arange(0, num_wlkrs.max(), dtype=float)
    shuf_orig = numpy.arange(num_wlkrs_max)
    rng.shuffle(shuf_orig)
    cshuf = struct()
    cshuf.ndata_last = None
    cshuf.order_orig = shuf_orig

  def coh_shuffle_data(ws1):
    """`Coherently' shuffle the dataset.
    """
    ndata_cur = len(ws1.X_w)
    if cshuf.ndata_last != ndata_cur:
      # Regenerates the shuffle data
      dbg(10, "coh_shuffle_data: regenerating order for ndata=%d\n" % ndata_cur)
      cshuf.order_new = shuffle_with_eliminations(ndata_cur, cshuf.order_orig)
      cshuf.ndata_last = ndata_cur
    order = cshuf.order_new
    X0, Xw0, w0 = ws1.X, ws1.X_w, ws1.w
    X1 = X0[order]
    Xw1 = Xw0[order]
    w1 = w0[order]
    ws.X = X1
    ws.X_w = Xw1
    ws.w = w1

  while True:
    raw_data = raw_dset.load_meas_data(stage=stage, match=match_proc)
    meas = raw_data.flat
    N_meas = len(meas)
    dbg(10, "%s %s %s %s %s %s\n" \
        % (N_meas, stage, match_proc, match_params['start'], match_params['stop'], i))

    # Prepare matcher for the next iteration
    match_params['start'] += match_params['skip']
    match_params['stop'] += match_params['skip']

    # Don't quit too early or else we will skip real data
    if N_meas == 0:
      continue

    if act_row_begin is None:
      # note the first and ending rows
      act_row_begin = meas[0].meta['row']
    act_row_end = meas[-1].meta['row']
    dbg(10, "act_row_begin, act_row_end = %d, %d\n" \
        % (act_row_begin, act_row_end))
    #if N_meas == 0:
    #  print "End of record: ", match_params
    #  break  # end-of-loop

    for m in meas:
      El_1 = m.El
      wt_1 = m.wtwlkr
      beta = m.meta['beta']
      # SEVERE DEBUGGING
      if __first:
        print wt_1
        __first = False
      # HACK PROCESSING (only if requested)
      if free_proj == "wrong-old-phaseless":
        wt = numpy.real(m.wtwlkr)
        El = numpy.real(m.El)
      # NORMAL PROCESSING BELOW
      elif free_proj:
        # "muck" these things for free projection
        wt = numpy.real(m.wtwlkr)
        El = numpy.real(m.wtwlkr * m.El) / wt
      else:
        wt = numpy.abs(m.wtwlkr)
        El = numpy.real(m.El)

      beta = m.meta['beta'] # imaginary time value
      Emix = numpy.average(El, weights=wt)
      #Emix_no_wt = numpy.average(El)
      #Emix_err = numpy.std(El) / numpy.sqrt(len(El) - 1)
      Emix_aa_jk = jk_generate_averages(a=El, weights=wt)
      (Emix_jk, Emix_jk_err, Emix_jk_corr) \
          = jk_stats_aa(aa_jk=Emix_aa_jk)
      # beware of order
      raw_meas[i] = (m.meta['count'], m.meta['beta'],
                     Emix, sum(wt), Emix_jk, Emix_jk_err)

      ws.init(El, wt)
      nsamp = m.meta['count']
      if pop_pre == 'shuffle':
        ws.shuffle(rng)
      elif pop_pre == 'coherent-shuffle':
        coh_shuffle_data(ws)
      for (pblk_count, pblk_size) in zip(pblkcounts, nsamp // pblkcounts):
        if pop_pre in ('none', 'contiguous', 'shuffle', 'coherent-shuffle'):
          ws.reblock_p(pblk_count, save=True)
        elif pop_pre in ('buggy-contiguous',):
          # WARNING: Can trigger bug if pop size is not commensurate
          # to pblk_count
          ws.reblock(pblk_size, save=True)
        elif pop_pre in ('buggy-transpose'):
          # WARNING: Can trigger bug if pop size is not commensurate
          # to pblk_count
          ws.reblock(pblk_size, transpose=True, save=True)
        #if pop_pre in ('none', 'contiguous', 'shuffle'):
        #  ws.reblock(pblk_size, save=True)
        #elif pop_pre in ('interleaving', 'transpose'):
        #  ws.reblock(pblk_size, transpose=True, save=True)
        else:
          raise ValueError, "Invalid pop_pre=%s" % str(pop_pre)
        rblk = ws._reblk_rec
        rblk_row = raw_planks[pblk_count][:,i]
        rblk_row['w'] = rblk.reblk_w
        rblk_row['Xw'] = rblk.reblk_Xw

      i = i + 1

    if act_row_end > match_row_end:
      raise RuntimeError, "Should not happen!"

    if act_row_end == match_row_end:
      # Marking the end of the loop
      dbg(10, "end of loop detected: %s\n" % act_row_end)
      break

  # Final: reblocking cross time axis
  for pblk_count in pblkcounts:
    for ip in xrange(pblk_count):
      d = raw_planks[pblk_count][ip,:]
      if numpy.any(d['w'] == 0): raise RuntimeError
      ws.init(d['Xw'] / d['w'], d['w'])
      for tblk_size in tblksizes:
        ws.reblock(tblk_size, save=True)
        raw_blocks[(pblk_count, tblk_size)][ip,:]['w'] = ws._reblk_rec.reblk_w
        raw_blocks[(pblk_count, tblk_size)][ip,:]['Xw'] = ws._reblk_rec.reblk_Xw

  #Fsum.close()

  return rslt


def h5meas_dump_reblk_raw_data(blk):
  """Dumps the reblocked datasets in a text table format.
  Used for debugging.

  Input `blk' object is the output of h5meas_reblk routine.
  """
  from wpylib.iofmt.text_output import text_output
  rawblk = blk.raw_blocks
  for (pblk, tblk) in rawblk.keys():
    fn = "dump-%dx%d" % (pblk, tblk)
    fd = text_output(fn)
    rblk = rawblk[pblk,tblk]
    X = rblk['Xw'] / rblk['w']
    w = rblk['w']
    for t in xrange(rblk.shape[1]):
      for p in xrange(rblk.shape[0]):
        fd("%3d %3d %16.9f %16.9f\n" % (p, t, X[p,t], w[p,t]))

    fd.close()


def h5meas_dump_reblk_stats(blk, fn):
  """Dumps the reblocking statistics in a text table format.
  Input `blk' object is the output of h5meas_reblk routine.
  The output is like this:

      #tblk   pblk           mean_g            err_g            mean_b            err_b              w_tot tbsz            mean_jk           err_jk
         25      1  -8092.438212090      0.003446285   -8092.437398380      0.003932366    77093.632027542    2    -8092.438213362      0.003502350
         10      1  -8092.438212090      0.003178238   -8092.437240991      0.003409631    77093.632027542    5    -8092.438220600      0.003732636
          5      1  -8092.438212090      0.003726210   -8092.437043970      0.004222428    77093.632027542   10    -8092.438238199      0.004836244
     ...

  The fields are:

    tblk: number of blocks in the imaginary-time axis

    tbsz: list of the number of blocks in the imaginary-time axis

    pblk: size of each block in the imaginary-time axis (or,
        the number of time slices averaged for reblocking along time axis)

    mean_g, err_g: the 'grand' weighted statistics mean and (biased) error estimate

    mean_b, err_b: the unweighted statistics mean and error estimate
        (sometimes called 'block average'--i.e. plain averaging
        of the reblocked data),

    mean_jk, err_jk: mean and error restimate from the delete-1 jackknife procedure.

    w_tot: total weight of all the raw data points accounted
        in the reblocking process.
  """

  from wpylib.iofmt.text_output import text_output
  from wpylib.text_tools import str_fmt_heading
  from wpylib.math.stats.jackknife1 import jk_generate_averages, jk_stats_aa
  from numpy import average, product, std, sqrt, sum
  rawblk = blk.raw_blocks
  #fn = "reblock-stats"
  fd = text_output(fn)
  fmt = "%5d %6d %16.9f %16.9f  %16.9f %16.9f   %16.9f %4d   %16.9f %16.9f\n"
  fd(str_fmt_heading(fmt) % ('#tblk', 'pblk', 'mean_g', 'err_g', 'mean_b', 'err_b', 'w_tot', 'tbsz', 'mean_jk', 'err_jk'))
  for (pblk, tblk) in sorted(rawblk.keys()):
    rblk = rawblk[pblk,tblk]
    ndata = product(rblk.shape)
    # reblocked "measurements":
    X_ds = rblk['Xw'] / rblk['w']
    wtot = sum(rblk['w'])
    w2tot = sum(rblk['w']**2)
    # grand (weighted) averages
    Xg = sum(rblk['Xw']) / sum(rblk['w'])
    # REF: http://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance
    # biased sample variance
    Xg_var_biased = sum(rblk['w'] * (X_ds - Xg)**2) / wtot
    #Xg_var_unbiased =
    # biased error estimate of the "grand" average
    Xg_err_biased = sqrt(Xg_var_biased / (ndata))
    # FIXME: Add also unbiased error estimate
    # poor but simplest attempt to reduce bias in var estimate above:
    Xg_err_unbiased1 = sqrt(Xg_var_biased / (ndata-1))
    # unweighted statistics: block average, standard error of the blk avg
    Xb = average(X_ds)
    Xerr = std(X_ds, ddof=1) / sqrt(ndata)
    # Jackknife statistics
    Emix_aa_jk = jk_generate_averages(a=X_ds.flatten(), weights=rblk['w'].flatten())
    (Emix_jk, Emix_jk_err, Emix_jk_corrected) = jk_stats_aa(aa_jk=Emix_aa_jk)

    fd(fmt % (rblk.shape[1], pblk, Xg, Xg_err_biased, \
              Xb, Xerr, wtot, tblk, Emix_jk, Emix_jk_err))
  fd.close()
