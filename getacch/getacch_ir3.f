c***********************************************************************
c                        GETACCH_IR3.F
c***********************************************************************
c Calculate r^-3 for an array of particles
c
c Input:
c    nbod     ==>  number of massive bodies
c    istart   ==>  body to start with
c    x,y,z    ==>  positions
c
c Output:
c    ir3      ==>  r^-3
c    ir       ==>  r^-1

      subroutine getacch_ir3(nbod, istart, x, y, z, ir3, ir)

      include '../swift.inc'

c...  Inputs:
      integer nbod, istart
      real*8 x(nbod), y(nbod), z(nbod)

c...  Outputs:
      real*8 ir3(nbod)
      real*8 ir(nbod)

c...  Internals:
      integer i
      real*8 r2

c----
c...  Executable code

      do i=istart,nbod
         r2 = x(i)*x(i) + y(i)*y(i) + z(i)*z(i)
         ir(i) = 1.0d0/sqrt(r2)
         ir3(i) = ir(i)/r2
      enddo

      return
      end       ! getacch_ir3
c-----------------------------------------------------------------------
