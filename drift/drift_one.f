c***********************************************************************
c                        DRIFT_ONE.F
c***********************************************************************
c This subroutine does the Danby-type drift for one particle, using
c appropriate variables and redoing a drift if the accuracy is too poor
c (as flagged by the integer iflg).
c
c Input:
c    nbod          ==>  Number of massive bodies
c    mu            ==>  Mass of central body
c    x, y, z       ==>  Initial position in Jacobi coord
c    vx, vy, vz    ==>  Initial velocity in Jacobi coord
c    dt            ==>  Time step
c
c Outputs:
c    x, y, z       ==>  Final position in Jacobi coord
c    vx, vy, vz    ==>  Final velocity in Jacobi coord
c    iflg          ==>  integer (zero for successful step)

      subroutine drift_one(mu, x, y, z, vx, vy, vz, dt, iflg)

      include '../swift.inc'

c...  Inputs Only:
      real*8 mu, dt

c...  Inputs and Outputs:
      real*8 x, y, z
      real*8 vx, vy, vz

c...  Output
      integer iflg

c...  Internals:
      integer i
      real*8 dttmp

c----
c...  Executable code

      call drift_dan(mu, x, y, z, vx, vy, vz, dt, iflg)

      if(iflg .ne. 0) then

         do i = 1,10
            dttmp = dt/10.d0
            call drift_dan(mu, x, y, z, vx, vy, vz, dttmp, iflg)
            if(iflg .ne. 0) return
         enddo

      endif

      return
      end    ! drift_one
c-----------------------------------------------------------------------
