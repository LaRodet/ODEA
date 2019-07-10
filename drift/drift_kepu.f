c***********************************************************************
c                        DRIFT_KEPU.F
c***********************************************************************
c Subroutine for solving kepler's equation using universal variables.
c
c Inputs:
c    dt            ==>  Time step
c    r0            ==>  Relative distance
c    mu            ==>  Reduced mass of system
c    alpha         ==>  Energy
c    u             ==>  Angular momentum
c
c Outputs:
c    fp            ==>  f' from p170
c    c1, c2, c3    ==>  c's from p171-172
c    iflg          ==>  =0 if converged; !=0 if not

      subroutine drift_kepu(dt, r0, mu, alpha, u, fp, c1, c2, c3, iflg)

      include '../swift.inc'

c...  Inputs:
      real*8 dt, r0, mu, alpha, u

c...  Outputs:
      real*8 fp, c1, c2, c3
      integer iflg

c...  Internals:
      real*8 s, st, fo, fn

c----
c...  Executable code

      call drift_kepu_guess(dt, r0, mu, alpha, u, s)

      st = s
c..   Store initial guess for possible use later in
c..   Laguerre's method, in case Newton's method fails.

      call drift_kepu_new(s, dt, r0, mu, alpha, u, fp, c1, c2, c3, iflg)
      if(iflg.ne.0) then
         call drift_kepu_fchk(dt, r0, mu, alpha, u, st, fo)
         call drift_kepu_fchk(dt, r0, mu, alpha, u, s, fn)
         if(abs(fo).lt.abs(fn)) then
            s = st
         endif
         call drift_kepu_lag(s, dt, r0, mu, alpha, u, fp, c1, c2, c3,
     &      iflg)
      endif

      return
      end    ! drift_kepu
c-----------------------------------------------------------------------
