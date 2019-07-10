c***********************************************************************
c                    ORBEL_EGET.F
c***********************************************************************
c PURPOSE:  Solves Kepler's equation.
c  For e < 0.18 use ESOLMD for speed and sufficient accuracy
c  For e > 0.8 use EHIE - this one may not converge fast enough.
c
c Inputs:
c    e ==> eccentricity anomaly
c    m ==> mean anomaly
c
c Returns:
c    orbel_eget ==> eccentric anomaly
c
c ALGORITHM: Quartic convergence from Danby
c REMARKS: For results very near roundoff, give it M between
c           0 and 2*pi. One can condition M before calling EGET
c           by calling my double precision function MOD2PI(M).
c           This is not done within the routine to speed it up
c           and because it works fine even for large M.

      real*8 function orbel_eget(e, m)

      include '../swift.inc'

c...  Inputs Only:
      real*8 e, m

c...  Internals:
      real*8 x, sm, cm, sx, cx
      real*8 es, ec, f, fp, fpp, fppp, dx

c----
c...  Executable code

      call orbel_scget(m, sm, cm)

c...  Begin with a guess accurate to order ecc**3
      x = m + e*sm*( 1.d0 + e*( cm + e*( 1.d0 -1.5d0*sm*sm)))

c... Go through one iteration for improved estimate
      call orbel_scget(x, sx, cx)
      es = e*sx
      ec = e*cx
      f = x - es  - m
      fp = 1.d0 - ec
      fpp = es
      fppp = ec
      dx = -f/fp
      dx = -f/(fp + dx*fpp/2.d0)
      dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)
      orbel_eget = x + dx

c... Do another iteration.
c... For m between 0 and 2*pi this seems to be enough to
c... get near roundoff error for eccentricities between 0 and 0.8
      x = orbel_eget
      call orbel_scget(x,sx,cx)
      es = e*sx
      ec = e*cx
      f = x - es  - m
      fp = 1.d0 - ec
      fpp = es
      fppp = ec
      dx = -f/fp
      dx = -f/(fp + dx*fpp/2.d0)
      dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)

      orbel_eget = x + dx

      return
      end  ! orbel_eget
c-----------------------------------------------------------------------
