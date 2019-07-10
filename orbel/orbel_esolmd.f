c***********************************************************************
c                    ORBEL_ESOLMD.F
c***********************************************************************
c PURPOSE:  Solves Kepler's equation.
c
c Inputs:
c    e ==> eccentricity anomaly
c    m ==> mean anomaly
c
c Returns:
c    orbel_esolmd ==> eccentric anomaly
c
c ALGORITHM: Some sort of quartic convergence from Wisdom.
c REMARKS: Only good for small eccentricity since it only
c         iterates once. (good for planet calcs.)
c         also does not put m or e between 0. and 2*pi

      real*8 function orbel_esolmd(e, m)

      include '../swift.inc'

c...  Inputs Only:
      real*8 e, m

c...  Internals:
      real*8 x, sm, cm, sx, cx
      real*8 es, ec, f, fp, fpp, fppp, dx

c----
c...  Executable code

      call orbel_scget(m, sm, cm)
      x = m + e*sm*( 1.d0 + e*( cm + e*( 1.d0 -1.5d0*sm*sm)))

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

      orbel_esolmd = x + dx

      return   ! orbel_esolmd
      end
c-----------------------------------------------------------------------
