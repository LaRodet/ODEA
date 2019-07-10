c***********************************************************************
c                        DRIFT_KEPU_P3SOLVE.F
c***********************************************************************
c Returns the real root of cubic often found in solving kepler
c  problem in universal variables.
c
c Inputs:
c    dt    ==>  time step (real scalar)
c    r0    ==>  Relative distance
c    mu    ==>  Reduced mass of system
c    alpha ==>  Twice the binding energy
c    u     ==>  Vel. dot radial vector
c
c Outputs:
c    s     ==>  solution of cubic eqn for the universal variable
c    iflg  ==>  success flag ( = 0 if O.K.)

      subroutine drift_kepu_p3solve(dt, r0, mu, alpha, u, s, iflg)

c...  Inputs:
      real*8 dt, r0, mu, alpha, u

c...  Outputs:
      integer iflg
      real*8 s

c...  Internals:
      real*8 denom, a0, a1, a2, q, r, sq2, sq, p1, p2

c----
c...  Executable code

      denom = (mu - alpha*r0)/6.d0
      a2 = 0.5*u/denom
      a1 = r0/denom
      a0 =-dt/denom

      q = (a1 - a2*a2/3.d0)/3.d0
      r = (a1*a2 -3.d0*a0)/6.d0 - (a2**3)/27.d0
      sq2 = q**3 + r**2

      if( sq2 .ge. 0.d0) then
         sq = sqrt(sq2)

         if ((r+sq) .le. 0.d0) then
            p1 =  -(-(r + sq))**(1.d0/3.d0)
         else
            p1 = (r + sq)**(1.d0/3.d0)
         endif
         if ((r-sq) .le. 0.d0) then
            p2 =  -(-(r - sq))**(1.d0/3.d0)
         else
            p2 = (r - sq)**(1.d0/3.d0)
         endif

         iflg = 0
         s = p1 + p2 - a2/3.d0

      else
         iflg = 1
         s = 0
      endif

      return
      end     !   drift_kepu_p3solve
c-----------------------------------------------------------------------
