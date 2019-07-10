c***********************************************************************
c                          ORBEL_ZGET.F
c***********************************************************************
c Purpose:  Solves the equivalent of Kepler's eqn. for a parabola
c  given Q (Fitz. notation.)
c
c Input:
c    q ==> parabola mean anomaly
c Returns:
c    orbel_zget ==> eccentric anomaly
c
c Algorithm: p. 70-72 of Fitzpatrick's book "Princ. of Cel. Mech."
c Remarks: For a parabola we can solve analytically.

      real*8 function orbel_zget(q)

      include '../swift.inc'

c...  Inputs Only:
      real*8 q

c...  Internals:
      integer iflag
      real*8 x, tmp

c----
c...  Executable code

      iflag = 0
      if(q.lt.0.d0) then
         iflag = 1
         q = -q
      endif

      if (q.lt.1.d-3) then
         orbel_zget = q*(1.d0 - (q*q/3.d0)*(1.d0 -q*q))
      else
         x = 0.5d0*(3.d0*q + sqrt(9.d0*(q**2) +4.d0))
         tmp = x**(1.d0/3.d0)
         orbel_zget = tmp - 1.d0/tmp
      endif

      if(iflag .eq.1) then
         orbel_zget = -orbel_zget
         q = -q
      endif

      return
      end    ! orbel_zget
c----------------------------------------------------------------------
