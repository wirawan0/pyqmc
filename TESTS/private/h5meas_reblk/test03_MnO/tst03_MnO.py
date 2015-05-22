# Imported 20150522 from Cr2_analysis_cbs.py
# (dated 20141017, CVS rev 1.143).

import numpy
from pyqmc.analysis.meas_reblk import h5_energy_meas_reblk as h5meas_reblk

#### TEST MEASUREMENT REBLOCKING SET 3: MnO ('Hybrid' Phaseless)

def test_mea_reblk3(show_summary=False, minbeta=0.0, maxbeta=9.5001, hack_reblk=False,
                    tblksizes=None, pblkcounts=None):
  """[20140310]

  Testcase:
    MnO/AFM2/rh.1x1x1/Opium-GFRG/vol11.14/dt005/k-0276-1744-0556.run/measurements.h5
    

  """
  from numpy import average, std, sqrt, sum
  from pyqmc.results.pwqmc_meas import meas_hdf5 as pwmeas_hdf5
  global modb3, moblk3

  if not "modb3" in globals():
    modb3 = pwmeas_hdf5("/home/wirawan/Work/PWQMC-77/expt/qmc/MnO/AFM2/rh.1x1x1/Opium-GFRG/vol11.14/dt005/k-0276-1744-0556.run/measurements.h5","r")
    # File fingerprint:
    # -rw-r----- 1 wirawan wirawan 14014256 2009-02-16 13:46:04 measurements.h5
    # f1887f95328f959788a8fb8736585519  measurements.h5
    # 8e7ceb2878c7500f9de7da2cb1550c88f1a87bdf  measurements.h5

  # opts
  pass
  # key values
  nblkstep = 100
  itv_Em = 20
  nblk_Em = nblkstep // itv_Em

  INFO_result = """
     0.5 5929.192791162419 -475.8849877487143
     1.0 5918.152949638847 -475.8916027288257
     1.5 5913.958045761421 -475.9818432947307
     2.0 5932.371149079349 -476.0211029719388
     2.5 5829.612095815706 -476.0058744868288
     3.0 5928.494407625599 -476.0598602392915
     3.5 5771.512233773923 -475.9219650550428
     4.0 5924.303884674163 -475.9059723036343
     4.5 5839.098729434451 -475.9055510204780
     5.0 5828.903657531735 -475.9340413602829
     5.5 5784.662578115800 -475.7920425875609
     6.0 5829.067921881763 -475.8217168957981
     6.5 5825.071917556909 -475.8353810046410
     7.0 5860.953482644119 -475.8748698192285
     7.5 5776.280049303745 -475.8742356220355
     8.0 5832.637528855600 -475.9218307507870
     8.5 5791.498419503304 -475.9859779863332
     9.0 5813.607630297450 -475.9458351372380
     9.5 5837.630550704131 -476.0496981979251
  """

  # Find out of internal reblocking will get us the same.
  #minbeta = 0
  #maxbeta = 9.5001
  def matchme(m):
    return (minbeta == None or m['beta'] >= minbeta) \
       and (maxbeta == None or m['beta'] <= maxbeta)

  if hack_reblk == "wrong-old-phaseless":
    free_proj = "wrong-old-phaseless"
  elif hack_reblk:
    raise NotImplementedError, "Not implemented: hack_reblk=%s" % (hack_reblk,)
  else:
    free_proj = False

  if not tblksizes:
    tblksizes = [nblk_Em]

  if not pblkcounts:
    pblkcounts = [1]

  moblk3 = h5meas_reblk(modb3, free_proj=free_proj,
               tblksizes=tblksizes,
               pblkcounts=pblkcounts,
               # Filter for records:
               match=matchme, stage=2,
               )

  if show_summary:
    print "Raw data:"
    for r in moblk3.raw_meas:
      beta, E, w = r[1], r[2], r[3]
      print beta, w, E

    pblk = pblkcounts[0]
    tblk = tblksizes[0]
    rawblk = moblk3.raw_blocks[pblk, tblk]
    #rawblk_ndata = rawblk.shape[1] # num of data across time axis
    rawblk_ndata = numpy.product(rawblk.shape) # num of data across ALL axes
    print "Reblocked data:"
    for (r,Xw,w) in zip(moblk3.raw_meas[tblk-1::tblk], rawblk['Xw'].T, rawblk['w'].T):
      print "%3.1f %.12f %.13f" % (r['beta'], sum(w), sum(Xw)/sum(w))

    #for (r,Xw,w) in zip(moblk3.raw_meas[tblk-1::tblk], rawblk['Xw'][0], rawblk['w'][0]):
    #  print "%3.1f %.12f %.13f" % (r['beta'], w, Xw/w)

    print "Grand average:", sum(rawblk['Xw']) / sum(rawblk['w'])
    print "Block average:", average(rawblk['Xw'] / rawblk['w'])
    print "Block error:  ", std(rawblk['Xw'] / rawblk['w'], ddof=1) / sqrt(rawblk_ndata)
    ANS = """
      In [567]: test_mea_reblk3(maxbeta=33.5001)
      matching: route 4
      Matching rows: 335 (30 .. 364)
      335 2 <function match_rows at 0x7f8e1410f500> 30 365 0
      act_row_begin, act_row_end = 30, 364
      end of loop detected: 364
      Reblocked data:
      0.5 5929.192791162424 -475.8849877487137
      1.0 5918.152949638845 -475.8916027288261
      1.5 5913.958045761422 -475.9818432947310
      2.0 5932.371149079355 -476.0211029719388
      2.5 5829.612095815713 -476.0058744868288
      3.0 5928.494407625602 -476.0598602392911
      ...
      Grand average: -475.904013131
      Block average: -475.903821338     <<
      Block error:   0.00939713762067   <<

  GOOD. These are reproducing analysis.groff in the meas dir above:

      ~/Work/PWQMC-77/expt/qmc/MnO/AFM2/rh.1x1x1/Opium-GFRG/vol11.14/dt005/k-0276-1744-0556.run $ check_reblocking raw.out  plt txt 1 20 1
      ...
           67       1      67      -475.9038213    0.005916515215    0.009397137715 -475.9038(94)
  """


def test_mea_reblk3b():
  """[20140310]
  Continuation of test_mea_reblk3().
  """
  global modb3, moblk3b
  moblk3b = h5meas_reblk(modb3, free_proj=False,
               tblksizes=[10],
               pblkcounts=[10],
               # Filter for records:
               match=None, stage=2,
               )


"""
In [546]: execfile("/home/wirawan/Work/GAFQMC/_scrap/scrap-20140310.py")

In [547]: test_mea_reblk3b()
matching: route 4
Matching rows: 338 (30 .. 367)
338 2 <function match_rows at 0xc9d2050> 30 368 0
act_row_begin, act_row_end = 30, 367
end of loop detected: 367

In [548]: X = sum( moblk3b.raw_blocks[10,11]['Xw'] ) / sum( moblk3b.raw_blocks[10,11]['w'])

In [549]: X
Out[549]: -475.90465671818532

In [551]: numpy.average( moblk3b.raw_blocks[10,11]['Xw'] / moblk3b.raw_blocks[10,11]['w'])
Out[551]: -475.90406123995518



In [614]: test_mea_reblk3(maxbeta=33.5001, pblkcounts=[10], tblksizes=[1])
matching: route 4
Matching rows: 335 (30 .. 364)
335 2 <function match_rows at 0x2833410> 30 365 0
act_row_begin, act_row_end = 30, 364
end of loop detected: 364
Reblocked data:
0.1 121.503598548362 -475.7795152937875
0.2 112.881478946709 -475.9734788765674
0.3 108.746858875414 -475.9194230019845
0.4 111.955027198162 -475.9720016211089
0.5 120.757204674677 -476.1074224178457
0.6 115.027654975497 -475.9572580101893
0.7 119.951165979957 -476.2674229638554
0.8 126.398664408384 -476.0150673799081
0.9 115.225913355990 -475.7918592559600
...
Grand average: -475.904114803
Block average: -475.901928561
Block error:   0.0144268116986
"""


def test_mea_reblk3c():
  """[20140311]
  Continuation of test_mea_reblk3().
  Should reproduce classic/legacy reblocking.
  """
  global modb3, moblk3c
  moblk3c = h5meas_reblk(modb3, free_proj=False,
               tblksizes=[5],
               pblkcounts=[1],
               # Filter for records:
               match=None, stage=2,
               )


def test_mea_reblk3d():
  """[20140311]
  Specialization of test_mea_reblk3().
  Attempt to reproduce

     /home/wirawan/Work/PWQMC-77/expt/qmc/MnO/AFM2/rh.1x1x1/Opium-GFRG/vol11.14/dt005/k-0276-1744-0556.run/reblocking-errbar.txt

  (no EL cap) to within the best we can reproduce it.
  Result: good, we can!
  """
  global modb3, moblk3, moblk3d

  #try:
  #  moblk3_save = moblk3
  #except:
  #  moblk3_save = None

  pblkcounts = numpy.array([100, 90, 80, 70, 60, 50, 40, 30, 20, 10])
  tblkcounts = numpy.array([67,33,16,8,4,2])
  # 338 is the row count of the loaded dataset
  tblksizes = 338 // tblkcounts
  print "pblkcounts = ", pblkcounts
  print "tblksizes  = ", tblksizes

  minbeta = 0.0
  maxbeta = 33.80001

  def matchme(m):
    return (minbeta == None or m['beta'] >= minbeta) \
       and (maxbeta == None or m['beta'] <= maxbeta)

  moblk3d = h5meas_reblk(modb3, free_proj=False,
               tblksizes=tblksizes,
               pblkcounts=pblkcounts,
               # Filter for records:
               match=matchme, stage=2,
               )

  #test_mea_reblk3(show_summary=False, minbeta=0.0, maxbeta=33.8001, hack_reblk=False,
  #                tblksizes=tblksizes, pblkcounts=pblkcounts)
  #moblk3d = moblk3

  #if moblk3_save != None:
  #  moblk3 = moblk3_save
  #else:
  #  del moblk3

  h5meas_dump_reblk_stats(moblk3d, "summary-reblocking.reblk3d.txt")
