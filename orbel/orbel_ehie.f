c***********************************************************************
c                    ORBEL_EHIE.F
c***********************************************************************
c PURPOSE:  Solves Kepler's equation.
c
c Inputs:
c    e ==> eccentricity anomaly
c    m ==> mean anomaly
c
c Returns:
c    orbel_ehybrid ==>  eccentric anomaly.
c
c ALGORITHM: Use Danby's quartic for 3 iterations.
c             Eqn. is f(x) = x - e*sin(x+M). Note  that
c             E = x + M. First guess is very good for e near 1.
c             Need to first get M between 0. and PI and use
c             symmetry to return right answer if M between PI and 2PI
c REMARKS: Modifies M so that both E and M are in range (0,TWOPI)

      real*8 function orbel_ehie(e, m)

      include '../swift.inc'

c...  Inputs Only:
      real*8 e, m

c...  Internals:
      integer iflag, nper, niter, NMAX
      real*8 dx, x, sa, ca, esa, eca, f, fp

      parameter (NMAX = 3)

c----
c...  Executable code

c...  In this section, bring M into the range (0,TWOPI) and if
c...  the result is greater than PI, solve for (TWOPI - M).
      iflag = 0
      nper = m/TWOPI
      m = m - nper*TWOPI
      if (m .lt. 0.d0) m = m + TWOPI

      if (m.gt.PI) then
         m = TWOPI - m
         iflag = 1
      endif

c...  Make a first guess that works well for e near 1.
      x = (6.d0*m)**(1.d0/3.d0) - m
      niter =0

c... Iteration loop
      do niter =1,NMAX
         call orbel_scget(x + m,sa,ca)
         esa = e*sa
         eca = e*ca
         f = x - esa
         fp = 1.d0 -eca
         dx = -f/fp
         dx = -f/(fp + 0.5d0*dx*esa)
         dx = -f/(fp + 0.5d0*dx*(esa+0.3333333333333333d0*eca*dx))
         x = x + dx
      enddo

      orbel_ehie = m + x

      if (iflag.eq.1) then
         orbel_ehie = TWOPI - orbel_ehie
         m = TWOPI - m
      endif

      return
      end         !orbel_ehie
c-----------------------------------------------------------------------
