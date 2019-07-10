c***********************************************************************
c                        DRIFT_KEPU_NEW.F
c***********************************************************************
c Solve Kepler's equation in universal variables using NEWTON's method
c
c Inputs:
c    s             ==>  Inital value of universal variable
c    dt            ==>  Time step
c    r0            ==>  Relative distance
c    mu            ==>  Reduced mass of system
c    alpha         ==>  Energy
c    u             ==>  Angular momentum
c
c Outputs:
c    s             ==>  Final value of universal variable
c    fp            ==>  f' from p170
c    c1, c2, c3    ==>  c's from p171-172
c    iflgn         ==>  =0 if converged; !=0 if not

      subroutine drift_kepu_new(s, dt, r0, mu, alpha, u, fp, c1, c2, c3,
     &   iflgn)

      include '../swift.inc'

c...  Inputs:
      real*8 s, dt, r0, mu, alpha, u

c...  Outputs:
      real*8 fp, c1, c2, c3
      integer iflgn

c...  Internals:
      integer nc
      real*8 x, c0, ds
      real*8 f, fpp, fppp, fdt

c----
c...  Executable code

      do nc=0,6
         x = s*s*alpha
         call drift_kepu_stumpff(x, c0, c1, c2, c3)
         c1 = c1*s
         c2 = c2*s*s
         c3 = c3*s*s*s
         f = r0*c1 + u*c2 + mu*c3 - dt
         fp = r0*c0 + u*c1 + mu*c2
         fpp = (-r0*alpha + mu)*c1 + u*c0
         fppp = (- r0*alpha + mu)*c0 - u*alpha*c1
         ds = - f/fp
         ds = - f/(fp + ds*fpp/2.0d0)
         ds = -f/(fp + ds*fpp/2.0d0 + ds*ds*fppp/6.0d0)
         s = s + ds
         fdt = f/dt

c..      Quartic convergence
         if( fdt*fdt.lt.DANBYB*DANBYB) then
             iflgn = 0
             return
         endif
c...     Newton's method succeeded

        enddo

c..     Newton's method failed
        iflgn = 1
        return

        end  ! drift_kepu_new
c-----------------------------------------------------------------------
