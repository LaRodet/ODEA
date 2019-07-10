c***********************************************************************
c                    ORBEL_EHYBRID.F
c***********************************************************************
c PURPOSE:  Solves Kepler's equation.
c
c Inputs:
c    e ==> eccentricity anomaly
c    m ==> mean anomaly
c
c Returns:
c    orbel_ehybrid ==> eccentric anomaly
c
c ALGORITHM: For e < 0.18 uses fast routine ESOLMD
c            For larger e but less than 0.8, uses EGET
c            For e > 0.8 uses EHIE
c REMARKS: Only EHIE brings M and E into range (0,TWOPI)

      real*8 function orbel_ehybrid(e,m)

      include '../swift.inc'

c...  Inputs Only:
      real*8 e,m

c...  Internals:
      real*8 orbel_esolmd, orbel_eget, orbel_ehie

c----
c...  Executable code

      if(e .lt. 0.18d0) then
         orbel_ehybrid = orbel_esolmd(e,m)
      else
         if( e .le. 0.8d0) then
            orbel_ehybrid = orbel_eget(e,m)
         else
            orbel_ehybrid = orbel_ehie(e,m)
         endif
      endif

      return
      end     ! orbel_ehybrid
c-----------------------------------------------------------------------
