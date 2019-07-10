c***********************************************************************
c                    ORBEL_FGET.F
c***********************************************************************
c PURPOSE: Solves Kepler's equation for hyperbola using hybrid approach
c
c Inputs:
c    e ==> eccentricity anomaly
c    capn ==> hyperbola mean anomaly
c
c Returns:
c     orbel_fget ==> eccentric anomaly
c
c ALGORITHM: Based on pp. 70-72 of Fitzpatrick's book "Principles of
c           Cel. Mech. ".  Quartic convergence from Danby's book.

      real*8 function orbel_fget(e, capn)

      include '../swift.inc'

c...  Inputs Only:
      real*8 e, capn

c...  Internals:
      integer i, IMAX
      real*8 tmp, x, shx, chx
      real*8 esh, ech, f, fp, fpp, fppp, dx
      PARAMETER (IMAX = 10)

c----
c...  Executable code

c... Begin with a guess proposed by Danby
      if( capn .lt. 0.d0) then
         tmp = -2.d0*capn/e + 1.8d0
         x = -log(tmp)
      else
         tmp = +2.d0*capn/e + 1.8d0
         x = log( tmp)
      endif

      orbel_fget = x

      do i = 1,IMAX
         call orbel_schget(x,shx,chx)
         esh = e*shx
         ech = e*chx
         f = esh - x - capn
         fp = ech - 1.d0
         fpp = esh
         fppp = ech
         dx = -f/fp
         dx = -f/(fp + dx*fpp/2.d0)
         dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)
         orbel_fget = x + dx
c...     If we have converged here there's no point in going on
         if(abs(dx) .le. TINY) RETURN
         x = orbel_fget
      enddo

      write(6,*) 'fget: returning without complete convergence'
      return
      end   ! orbel_fget
c-----------------------------------------------------------------------
