c***********************************************************************
c                            CE_CRIT.F
c***********************************************************************
c This subroutine computes criteria to decide if a hierarchy change is needed
c Only does massive particles
c
c Inputs:
c    nbod          ==> Number of massive bodies
c    mass          ==> Masses of bodies
c    oloc          ==> Link Bodies <--> Orbits
c    eta, mu       ==> Masses for centers & sats
c    xb, yb, zb    ==> Position in Barycentric coord
c    xj, yj, zj    ==> Position in Gen. Jacobi coord
c    vxj, vyj, vzj ==> Velocity in Gen. Jacobi coord
c    axj, ayj, azj ==> Acceleration in Gen. Jacobi coord
c    dt            ==> Time step
c
c Output:
c    checkchange   ==> True if a hierarchy check is needed

      subroutine ce_crit(nbod, xj, yj, zj, mass, mu, eta, oloc,
     &     axj, ayj, azj, checkchange)

      include '../swift.inc'
      include '../odea.inc'

c...  Inputs Only:
      integer nbod, oloc(NPLMAX,NPLMAX)
      real*8 mass(nbod), eta(nbod), mu(nbod)
      real*8 xj(nbod), yj(nbod), zj(nbod)
      real*8 axj(nbod), ayj(nbod), azj(nbod)

c...  Outputs Only
      logical checkchange

c...  Internals:
      integer i, j, k, l
      real*8 r2j(nbod), akep, c(nbod)
      logical ok

c----
c...  Executable code

      do k=2,nbod
         r2j(k) = xj(k)*xj(k)+yj(k)*yj(k)+zj(k)*zj(k)
      end do

      checkchange = .False.
      do k=2, nbod
         akep = (mu(k)+eta(k))/r2j(k)
         c(k) = sqrt(axj(k)*axj(k)+ayj(k)*ayj(k)+azj(k)*azj(k))/akep
         if (c(k).gt.CRITCHANGE) checkchange = .True.
c     print*,k,c(k)
      end do

      return
      end   ! ce_crit
c-----------------------------------------------------------------------
