#!/usr/bin/python
# $Id: datfile.py,v 1.4 2011-09-07 14:26:00 wirawan Exp $
#
# pyqmc.gamess.datfile module
#
# Wirawan Purwanto
# Created: 20110901
#

import numpy
import os
import os.path
import sys

import wpylib.shell_tools as sh
import wpylib.regexps
from wpylib.sugar import ifelse
from wpylib.iofmt.text_output import text_output
from wpylib.iofmt.text_input import text_input

from pyqmc import PyqmcDataError


class movecs(object):
  """Basic representation of MO vector block reader for GAMESS *.dat file.
  """
  def init(self, *_argsl, **_argsd):
    if _argsl or _argsd:
      self.read(*argsl, **argsd)
  def read(self, infile, vec_select=1, verbose=0, nbasis=None, out=sys.stdout):
    """Reads off molecular orbital vectors.
    Usage:
      movecs = pyqmc.gamess.datfile.movecs(fname, [options])
    Valid options:
      vec_select = <integer>  (default: 1; 1-based choice)
      verbose    = 0|1        (default: 0)
      nbasis     = <integer>  (default: autodetected)
   
    This routine was translated from Gamess::ReadGamessMOVecs routine
    in my Gamess.pm perl module.
    The latter routine was derived from C2_UHF_gamess.pl dated ~20070813.
   
    CAUTION:
    The resulting orbital (orbitals_alpha, orbitals_beta) arrays are 1-based,
    both in the orbital index and the basis index.

    Strict vector ordering (1, 2, 3, ..., N) is required in the $VEC data.
    We will check orbital indices strictly.
    This requires the orbitals be strictly ordered, with no skipping, etc.
    Strict checking is necessary for proper reading when we have more than 99
    basis functions.
   
    In anticipating for large basis size, the rule for deducing UHF-type
    movecs is more complicated. For nbasis >= 100, the tag number rolls back
    to zero, unfortunately, which makes thing a bit difficult: when we see a
    tag of " 1" again, is it UHF beta sector, or movec #101?
    One way we can be assured that it IS an UHF movecs is prohibiting the size
    of movecs to be greater than the deduced nbasis, which is a reasonable
    restriction. Then, when we apparently encounter movec "101" when nbasis==100,
    we can be sure that the 101st vector is actually beta movec #1.
    Thus UHF movecs can be detected by the following rule:
    - ( old $VecTag != 0, or old $VecIndex == $nbasis already ) AND
      new $VecTag == 1 .
   
    FIXME:
    The solution above is still not foolproof in two cases:
   
    1) suppose we have an UHF-type movecs with nbasis=220,
    but each spin sector only has 100 orbitals listed.
    Then this will be interpreted as an RHF-type movecs with nbasis=220
    and norbitals=200.
    2) in spherical basis, maximum norbitals is <= nbasis.
    When this happens, then the deduced nbasis is not the right number of
    spherical basis functions (thus the nbasis deduced above is not right).
   
    But I haven't seen this case yet, so forget about them temporarily.
    """
    from wpylib.regexps import regex
    # MOVECS comments (always 3 lines preceding the $VEC block)
    rx_comment = regex(r'^--- ')
    rx_vec = regex(r'^(?i) \$vec')
    rx_endvec = regex(r'^(?i)  ?\$end')
    F = text_input(infile)
    comments = None
    found = False
    vec_blk_count = 0
    n_comment_lines = 0
    O = text_output(out, flush=True)
    spin = "alpha"
    udet = False
    AllVecs = {}
    for txt in F:
      if rx_comment % txt:
        comments = []
        n_comment_lines = 3
      if n_comment_lines > 0:
        comments.append(txt)
        n_comment_lines -= 1
      if rx_vec % txt:
        vec_blk_count += 1
        if vec_blk_count < vec_select:
          # the comments we just read (if any) are irrelevant, so
          # remove them.
          comments = None
          continue

        found = True
        # This is the actual movecs reading loop ---
        # The $END marker for initial orbital guess (PUNMO=.TRUE.) is
        # buggy--we must tolerate that
        txt = F.next()
        # VecIndex = MO index to identify the whole vector
        # VecTag = MO "tag" number
        # In general VecTag is equal to VecIndex except when we have >= 100
        # basis funcs (where VecTag has only the last two digits).
        # NOTE: VecTag is *always* a 2-character string!
        VecIndex = 1
        VecTag = " 1"
        AmplIndex = 0
        Ampl = []
        Vecs = [Ampl]
        while not (rx_endvec % txt):
          NewVecTag = txt[0:2]
          #print "H: $txt\n";
          #print "V: $NewVecTag\n";

          # We should safely assume that VecTag > 1 at the end
          # of alpha orbitals; if that's not the case, that's
          # YOUR fault (why doing 1-basis quantum chemistry?)
          if NewVecTag != VecTag:
            # Just in case, we are very pedantic in checking for errors here:
            if nbasis != None:
              # (1) nbasis must be consistent
              if AmplIndex != nbasis:
                raise PyqmcDataError, \
                  ("%s:%d: Inconsistent nbasis " + \
                   "(original guess was = %d, most recently deduced value = %d) " + \
                   "for %s vector #%d") \
                  % (infile, F.lineno, nbasis, AmplIndex, spin, VecIndex)
            else:
              # Deduce nbasis
              nbasis = AmplIndex
              if nbasis == 0:
                raise PyqmcDataError, \
                  ("%s:%d: nbasis detected as zero! " + \
                   "Maybe there is corruption in the input file?") \
                  % (infile, F.lineno)
              if verbose > 0:
                O("pyqmc.gamess.movecs.read:%s: Deduced nbasis = %d\n" \
                  % (infile, nbasis))

            # UHF-type vector detection scheme: see the notes above
            if (VecTag != " 0" or VecIndex == nbasis) and NewVecTag == " 1":
              if verbose > 0:
                O("pyqmc.gamess.movecs.read:%s: Found UHF-type movecs\n" \
                  % (infile,))

              if udet:
                raise PyqmcDataError, \
                  ("%s:%d: alpha and beta orbitals were already defined?! " + \
                   "Maybe there is a mistake with your $VEC data?") \
                  % (infile, F.lineno)

              AllVecs[spin] = numpy.array(Vecs, dtype=float).T
              # start all over with a new MO block
              Ampl = []
              Vecs = [Ampl]
              spin = "beta"
              udet = True
              VecIndex = 0
              # end UHF-type detection scheme
            else:
              # Some additional error checking(s)
              if VecIndex >= nbasis: # and nbasis != 100:
                # NOTE: We disallow norbitals > nbasis in our reader.
                raise PyqmcDataError, \
                  ("%s:%d: The $VEC block has more than %d orbitals, " + \
                   "which is prohibited by this routine") \
                  % (infile, F.lineno, nbasis)

              Ampl = []
              Vecs.append(Ampl)

            AmplIndex = 0; # Start off a new vector
            VecIndex += 1
            VecTag = NewVecTag
          # end new vector/new spin sector detection

          # Strict index vs. tag checking:
          TagChk = "%2d" % (VecIndex % 100)
          if TagChk != VecTag:
            raise PyqmcDataError, \
              ("%s:%d: Mismatch vector tag number in vector #%d " + \
               "(wanted '%s', got '%s')") \
              % (infile, F.lineno, VecIndex, TagChk, VecTag)

          # the amplitudes are always stored in -n.nnnnnnnnE+nn fixed format
          # (15 characters wide)
          txtdata = txt[5:].rstrip()
          lendata = len(txtdata) // 15

          Ampl += [ float(txtdata[i*15:i*15+15]) for i in xrange(0, lendata) ]
          AmplIndex += lendata

          # TODO: $i < 5 is allowed ONLY on the last line;
          # Make sure we catch that.
          #print $VecIndex, " ", $AmplIndex, "\n";

          try:
            txt = F.next()
          except StopIteration:
            raise PyqmcDataError, \
              ("%s:%d: Unexpected EOF while reading in $VEC data") \
              % (infile, F.lineno)
        # end loop for reading in $VEC data

        # Finalization: do final checks, etc.

        AllVecs[spin] = numpy.array(Vecs, dtype=float).T

        if AmplIndex != nbasis:
          raise PyqmcDataError, \
            ("%s:%d: Inconsistent nbasis " + \
             "(original guess was = %d, most recently deduced value = %d) " + \
             "for %s vector #%d") \
            % (infile, F.lineno, nbasis, AmplIndex, spin, VecIndex)

        if udet:
          if AllVecs['alpha'].shape != AllVecs['beta'].shape:
            raise PyqmcDataError, \
              ("%s:%d: Inconsistent shape of MO matrix: " + \
               "(alpha = %s, beta = %s)") \
              % (infile, F.lineno, \
                 AllVecs['alpha'].shape, \
                 AllVecs['beta'].shape, \
                )
        if verbose > 0:
          O("pyqmc.gamess.movecs.read:%s: Total MO vectors read = %s%s\n" \
            % (infile, VecIndex, ifelse(udet, " (UHF-type)", "")))

        # stop reading if the desired vectors have been loaded
        break

    # end main text reading loop

    if not found:
      raise PyqmcDataError, \
        ("%s: Cannot find $VEC block number %s") \
        % (infile, vec_select)

    # Save the reading results to the "self" record:
    self.filename = infile
    self.vec_select = vec_select
    self.comments = comments
    self.udet = udet
    self.nbasis = nbasis
    for (spin, mo) in AllVecs.iteritems():
      setattr(self, spin, mo)
    return self

  def str(self):
    rslt = [ " $VEC\n" ]
    for spin in ifelse(self.udet, ('alpha', 'beta'), ('alpha',)):
      vecs = getattr(self, spin)
      for (i,v) in enumerate(vecs.T):
        for j1 in xrange(0, self.nbasis, 5):
          j2 = min(j1+5, self.nbasis)
          rslt.append(("%2d%3d" % (i + 1, j1//5 + 1)) \
                      + "".join([ "%15.8E" % v[j] for j in xrange(j1,j2)]) \
                      + "\n")
    rslt.append(" $END\n")
    return "".join(rslt)

  def write(self, outfile):
    """Writes molecular orbital in GAMESS format.
    What written depends on the `nbasis` and `udet` attributes, and
    the data is in the `alpha` and (optionally) `beta` attributes.
    Comments are not written out; they must be appended manually if you
    want them.
    This method use the flexible `text_output` facility, so the outfile
    can be an open file object or a filename.
    """
    F = text_output(outfile)
    F(self.str())
    F.flush()


