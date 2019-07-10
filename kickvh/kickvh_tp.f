c***********************************************************************
c                        KICKVH_TP.F
c***********************************************************************
c To kick the velocity components v by a*dt for test particles
c
c Inputs:
c    ntp           ==>  Number of tps
c    vxt, vyt, vzt ==>  Initial velocity (Jacobi or heliocentric)
c    axt, ayt, azt ==>  acceleration (Jacobi or heliocentric)
c    istat         ==>  status of the test particles
c                      istat(i) = 0 ==> active:  = 1 not
c             NOTE: it is really a 2d array but we only use the 1st row
c                          dt   ==>  time step
c Outputs:
c    vxt, vyt, vzt ==>  Final velocity (Jacobi or heliocentric)

      subroutine kickvh_tp(ntp, vxt, vyt, vzt, axt, ayt, azt, istat, dt)

      include '../swift.inc'

c...  Inputs Only:
      integer ntp,istat(ntp)
      real*8 axt(ntp), ayt(ntp), azt(ntp)
      real*8 dt

c...   Inputs and Output:
      real*8 vxt(ntp), vyt(ntp), vzt(ntp)

c...  Internals:
      integer n

c----
c...  Executable code

      do n= 1, ntp
         if(istat(n).eq.0) then
            vxt(n) = vxt(n) + axt(n)*dt
            vyt(n) = vyt(n) + ayt(n)*dt
            vzt(n) = vzt(n) + azt(n)*dt
         endif
      enddo

      return
      end    ! kickvh_tp
c-----------------------------------------------------------------------
