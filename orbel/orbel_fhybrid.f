c***********************************************************************
c                    ORBEL_FHYBRID.F
c***********************************************************************
c PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.
c
c Inputs:
c    e ==> eccentricity anomaly
c    n ==> hyperbola mean anomaly
c
c Returns:
c    orbel_fhybrid ==> eccentric anomaly
c
c ALGORITHM: For abs(N) < 0.636*ecc -0.6 , use FLON
c            For larger N, uses FGET

      real*8 function orbel_fhybrid(e, n)

      include '../swift.inc'

c...  Inputs Only:
      real*8 e, n

c...  Internals:
      real*8 abn
      real*8 orbel_flon, orbel_fget

c----
c...  Executable code

      abn = n
      if(n.lt.0.d0) abn = -abn

      if(abn .lt. 0.636d0*e -0.6d0) then
         orbel_fhybrid = orbel_flon(e, n)
      else
         orbel_fhybrid = orbel_fget(e, n)
      endif

      return
      end  ! orbel_fhybrid
c-----------------------------------------------------------------------
