c***********************************************************************
c                        KICKVH.F
c***********************************************************************
c Kick the velocity components v by a*dt
c
c Inputs:
c    nbod         ==>  Number of bodies
c    vx, vy, vz   ==>  Initial velocity (Jacobi or heliocentric)
c    ax, ay, az   ==>  Acceleration (Jacobi or heliocentric)
c    dt           ==>  Time step
c
c Outputs:
c    vx, vy, vz   ==>  Final velocity (Jacobi or heliocentric)

      subroutine kickvh(nbod, vx, vy, vz, ax, ay, az, dt)

      include '../swift.inc'

c...  Inputs Only:
      integer nbod
      real*8 ax(nbod), ay(nbod), az(nbod)
      real*8 dt

c...  Inputs and Output:
      real*8 vx(nbod), vy(nbod), vz(nbod)

c...  Internals:
      integer n

c----
c...  Executable code

      do n= 2, nbod
         vx(n) = vx(n) + ax(n)*dt
         vy(n) = vy(n) + ay(n)*dt
         vz(n) = vz(n) + az(n)*dt
      enddo

      return
      end    ! kickvh
c-----------------------------------------------------------------------
