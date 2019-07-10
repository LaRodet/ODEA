c***********************************************************************
c                        DRIFT_TP.F
c***********************************************************************
c Compute the independent Keplerian motions of the tps
c  with the Danby routine
c
c Inputs:
c    ntp              ==>  Number of test particles
c    mcen             ==>  Mass of the centers
c    xjt, yjt, zjt    ==>  Initial position in Jacobi coord
c    vxjt, vyjt, vzjt ==>  Initial velocity in Jacobi coord
c    istat            ==>  Status of the test paricles
c                           istat(i,1) = 0 ==> active:  = 1 not
c                           istat(i,2) = -1 ==> Danby did not work
c    dt               ==>  Time step
c
c Outputs:
c    xjt, yjt, zjt    ==>  Final position in Jacobi coord
c    vxjt, vyjt, vzjt ==>  Final velocity in Jacobi coord

      subroutine drift_tp(ntp, mcen, xjt, yjt, zjt, vxjt, vyjt, vzjt,
     &   dt, istat)

      include '../swift.inc'

c...  Inputs Only:
      integer ntp
      real*8 mcen, dt

c...  Inputs and Outputs:
      integer istat(NTPMAX,NSTAT)
      real*8 xjt(ntp), yjt(ntp), zjt(ntp)
      real*8 vxjt(ntp), vyjt(ntp), vzjt(ntp)

c...  Internals:
      integer j,iflg

c----
c...  Executable code

c...  Take a drift forward dt
      do j = 1,ntp
         if(istat(j,1).eq.0) then
            call drift_one(mcen, xjt(j), yjt(j), zjt(j),
     &             vxjt(j), vyjt(j), vzjt(j), dt, iflg)
            if(iflg.ne.0) then
               istat(j,1) = 1
               istat(j,2) = -1
            endif
         endif
      enddo

      return
      end ! drift_tp
c-----------------------------------------------------------------------
