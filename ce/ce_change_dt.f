c***********************************************************************
c                            CE_CHANGE_DT.F
c***********************************************************************
c This subroutine updates the time step after hierarchy change
c Only does massive particles
c
c Inputs:
c    nbod          ==> Number of massive bodies
c    eta, mu       ==> Masses for centers & sats
c    xj, yj, zj    ==> Position in Gen. Jacobi coord (current hierarchy)
c    vxj, vyj, vzj ==> Velocity in Gen. Jacobi coord (current hierarchy)

c
c Outputs:
c   dt            ==> Time step

      subroutine ce_change_dt(nbod, eta, mu, xj, yj, zj, vxj, vyj, vzj,
     &   dt)

      include '../swift.inc'

c...  Inputs Only:
      integer nbod

c...  Inputs and Outputs:
      real*8 eta(nbod), mu(nbod)
      real*8 xj(nbod), yj(nbod), zj(nbod)
      real*8 vxj(nbod), vyj(nbod), vzj(nbod)

c...  Outputs Only
      real*8 dt

c...  Internals:
      integer k, ialpha
      real*8 sma, e, capm, T, inc, omega, capom, mtot

c----
c...  Executable code

      dt = 0.d0
      do k=2,nbod
         mtot = mu(k)+eta(k)
         call orbel_xv2el(xj(k), yj(k), zj(k), vxj(k), vyj(k), vzj(k),
     &      mtot, ialpha, sma, e, inc, capom, omega, capm)
         if (ialpha.eq.-1) then
            T = twopi*sqrt((sma*(1-e))**3/mtot)
c            T = twopi*sqrt(sma**3/mtot)
         else
            if (capm.lt.0.) then
               T = twopi*sqrt((sma*(e-1))**3/mtot)
            else
               T = twopi*sqrt(sqrt(xj(k)*xj(k)+yj(k)*yj(k)
     &              +zj(k)*zj(k))**3/mtot)
            end if
         end if
         if (k.eq.2) then
            dt = T/100.d0
         else
            dt = min(dt,T/100.d0)
         end if
      end do

      return
      end  ! ce_change_dt
c-----------------------------------------------------------------------
