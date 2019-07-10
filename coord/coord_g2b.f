c***********************************************************************
c	                    COORD_G2B.F
c***********************************************************************
c Converts from Generalized Jacobi to Barycentric coordinates.
c
c Inputs:
c    nbod ==> Number of bodies
c    umat ==> Passage matrix
c    mass ==> Masses of the bodies
c    xj, yj, zj ==> Generalized Jacobi coordinates
c    vxj, vyj, vzj ==> Generalized Jacobi velocities
c
c Outputs:
c    xb, yb, zb ==> Barycentric coordinates
c    vxb, vyb, vzb ==> Barycentric velocities

      subroutine coord_g2b(nbod, umat, mass, xj, yj, zj, vxj, vyj, vzj,
     &   xb, yb, zb, vxb, vyb, vzb)

      include '../swift.inc'

c...  Inputs:
      integer nbod
      real*8 mass(nbod), umat(NPLMAX,NPLMAX)
      real*8 xj(nbod), yj(nbod), zj(nbod)
      real*8 vxj(nbod), vyj(nbod), vzj(nbod)

c...  Outputs:
      real*8 xb(nbod), yb(nbod), zb(nbod)
      real*8 vxb(nbod), vyb(nbod), vzb(nbod)

c...  Internals:
      integer i, j

c----
c...  Executable code

c...  Matrix product      XB=UMAT*XJ
      do i = 1,nbod
         xb(i) = 0.0d0
         yb(i) = 0.0d0
         zb(i) = 0.0d0
         vxb(i) = 0.0d0
         vyb(i) = 0.0d0
         vzb(i) = 0.0d0
         do j = 2,nbod   ! 2 because first jac. coord is zero
            if (umat(i,j).ne.0.0d0) then
               xb(i) = xb(i) + umat(i,j)*xj(j)
               yb(i) = yb(i) + umat(i,j)*yj(j)
               zb(i) = zb(i) + umat(i,j)*zj(j)
               vxb(i) = vxb(i) + umat(i,j)*vxj(j)
               vyb(i) = vyb(i) + umat(i,j)*vyj(j)
               vzb(i) = vzb(i) + umat(i,j)*vzj(j)
            end if
         end do
      end do

      return
      end    ! coord_g2b
c-----------------------------------------------------------------------
