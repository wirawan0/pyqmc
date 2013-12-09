

class gafqmc_stat(object):
  '''GAFQMC statistics file.'''
  # Timing fields which are averaged over time:
  timing_fields = ('step', 'ovlp', 'CalcFB', 'ExpH', 'Elocal', 'gasdev')
  # stat_dtype0 is the structure of tables in the pwaf-*.stat files
  stat_dtype0 = [
    ('idx', int),
    ('T_elapsed', float),
  ]
  for nnn in timing_fields:
    stat_dtype0 += [ ('t_' + nnn, float), ('n_' + nnn, int) ]
  stat_dtype = numpy.dtype(stat_dtype0)


def gaf_stat_load(file, decumulate=1):
  '''Loads a single gaf-NNNNN.stat ("fort.17") file into an array.
  The output records are still cumulative in nature:
  No postprocessing other than laying out the data in structured form is done.'''
  F = open_input_file(file)
  rslt0 = []
  # Remove all NaNs and return them as zero
  NaN_to_zero = lambda x : x if not numpy.isnan(x) else 0
  for ll in F:
    L = ll.strip()
    # Right now simply use this criterion as a marker of a record, i.e.
    # the first field being a number: it's not foolproof, man:
    if L[0].isdigit():
      TV = zip(pwqmc_stat.stat_dtype0, L.split())
      rslt0.append( tuple([ NaN_to_zero(typx[1](val)) for (typx, val) in TV ]) )

  F.close()
  rslt = numpy.array(rslt0, dtype=gafqmc_stat.stat_dtype)
  if decumulate:
    return gaf_stat_decumulate(rslt)
  else:
    return rslt


def gaf_stat_decumulate(rec):
  '''Converts the stat record from the cumulative (original) to per-block quantities.'''
  rslt = rec.copy()
  #print rslt
  t_elapsed = rslt['T_elapsed']
  # Unfortunately numpy does not have a converse function of "accumulate".
  # Work on them in batch as much as possible:
  Nfields = len(pwqmc_stat.timing_fields)
  tval = numpy.empty((Nfields, len(rec)), dtype="=f8")
  nval = numpy.empty((Nfields, len(rec)), dtype=int)
  for (i,col) in zip(xrange(Nfields),pwqmc_stat.timing_fields):
    tval[i] = rec['t_' + col]
    nval[i] = rec['n_' + col]
  #print tval; print nval
  # make tval a sum instead of the accumulated average
  tval *= nval
  #print tval
  # "de"-cumulate the sums
  for i in xrange(len(rec)-1, 0, -1):
    tval[:,i] -= tval[:,i-1]
    nval[:,i] -= nval[:,i-1]
    t_elapsed[i] -= t_elapsed[i-1]
  nval_nz = nval.copy()
  #print tval; print nval
  # Avoid division by zero:
  numpy.putmask(nval_nz, nval_nz==0, 1)
  # reconstruct the block averages
  tval /= nval_nz
  #print nval_nz
  # Reformat into the correct array form:
  for (i,col) in zip(xrange(Nfields),gafqmc_stat.timing_fields):
    rslt['t_' + col] = tval[i]
    rslt['n_' + col] = nval[i]
  return rslt


