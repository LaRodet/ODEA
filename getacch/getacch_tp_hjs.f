c***********************************************************************
c                        GETACCH_TP_HJS.F
c***********************************************************************
c This subroutine calculates the acceleration on one TP in the HJS case
c
c Input:
c    nbod          ==>  number of massive bodies
c    mass          ==>  mass of bodies
c    eta           ==> Masses of centers for orbits
c    mu            ==> Masses of satellites for orbits
c    xj, yj, zj    ==> Pl. positions in jacobi coord
c    xb, yb, zb    ==> Pl. positions in bary coord
c    ir3j          ==> Inverse Jacobi radii^3
c    oloct         ==> Link between tp and orbits
c    etatp         ==> Masses of center for tp
c    xjt, yjt, zjt ==>  tp position in jacobi coord
c    xbt, ybt, zbt ==>  tp position in bary coord
c
c Output:
c    axbt, aybt, azbt ==> tp accel. in bary. coord (real arrays)

      subroutine getacch_tp_hjs(nbod, mass, eta, mu, xj, yj, zj,
     &     xb, yb, zb, ir3j, oloct, etatp, xbt, ybt, zbt, xjt, yjt, zjt,
     &     axbt, aybt, azbt)

      include '../swift.inc'

c...  Inputs:
      integer nbod
      integer oloct(nbod)
      real*8 mass(nbod), eta(nbod), mu(nbod)
      real*8 xj(nbod), yj(nbod), zj(nbod), ir3j(nbod)
      real*8 xb(nbod), yb(nbod), zb(nbod)
      real*8 xjt, yjt, zjt, xbt, ybt, zbt, etatp

c...  Outputs:
      real*8 axbt, aybt, azbt

c...  Internals:
      integer i, j
      real*8 fac, dx, dy, dz
      real*8 rji2, ir3jt, irij3, irjt

c----
c...  Executable code

c...  get the r^-3
      call getacch_ir3(1, 1, xjt, yjt, zjt, ir3jt, irjt)

c... The first term relative to the orbit of the tp
      fac = etatp*ir3jt
      axbt = fac*xjt
      aybt = fac*yjt
      azbt = fac*zjt

c...  now the jacobi terms
      do j=2,nbod    !  Check all the orbits
        if (oloct(j).eq.-1) then !  tp is a center in orbit #i
          fac = -mu(j)*ir3j(j)
          axbt = axbt + fac*xj(j)
          aybt = aybt + fac*yj(j)
          azbt = azbt + fac*zj(j)
        else if (oloct(j).eq.1) then ! tp is a satellite in orbit #i
          fac = eta(j)*ir3j(j)
          axbt = axbt + fac*xj(j)
          aybt = aybt + fac*yj(j)
          azbt = azbt + fac*zj(j)
        end if
      end do

c...  now the third terms. We need to consider the bodies

      do i=1,nbod
        dx = xbt - xb(i)
        dy = ybt - yb(i)
        dz = zbt - zb(i)
        rji2 = dx*dx + dy*dy + dz*dz
        irij3 = 1.0d0/(rji2*sqrt(rji2))
        fac = mass(i)*irij3

        axbt = axbt - fac*dx
        aybt = aybt - fac*dy
        azbt = azbt - fac*dz
      end do

      return
      end      ! getacch_tp_hjs
c-----------------------------------------------------------------------
