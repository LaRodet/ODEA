c***********************************************************************
c                        DRIFT_KEPU_GUESS.F
c***********************************************************************
c Initial guess for solving kepler's equation using universal variables.
c
c Inputs:
c    dt    ==>  Time step
c    r0    ==>  Relative distance
c    mu    ==>  Reduced mass of system
c    alpha ==>  Energy
c    u     ==>  Angular momentum
c
c Output:
c    s     ==>  Initial guess for the value of universal variable

      subroutine drift_kepu_guess(dt, r0, mu, alpha, u, s)

      include '../swift.inc'

c...  Inputs:
      real*8 dt, r0, mu, alpha, u

c...  Output:
      real*8 s

c...  Internals:
      integer iflg
      real*8 y, sy, cy, sigma, es
      real*8 x, a
      real*8 en, ec, e

c----
c...  Executable code

      if (alpha.gt.0.0) then
c...     Find initial guess for elliptic motion
         if( dt/r0 .le. 0.4)  then
            s = dt/r0 - (dt*dt*u)/(2.0*r0*r0*r0)
            return
         else
            a = mu/alpha
            en = sqrt(mu/(a*a*a))
            ec = 1.0 - r0/a
            es = u/(en*a*a)
            e = sqrt(ec*ec + es*es)
            y = en*dt - es
            call orbel_scget(y, sy, cy)
            sigma = dsign(1.d0,(es*cy + ec*sy))
            x = y + sigma*.85*e
            s = x/sqrt(alpha)
         endif

      else
c...     Find initial guess for hyperbolic motion.
         call drift_kepu_p3solve(dt, r0, mu, alpha, u, s, iflg)
         if(iflg.ne.0) then
            s = dt/r0
         endif
      endif

      return
      end     !   drift_kepu_guess
c-----------------------------------------------------------------------
