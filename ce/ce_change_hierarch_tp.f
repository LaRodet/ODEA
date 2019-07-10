c***********************************************************************
c                            CE_CHANGE_HIERARCH_TP.F
c***********************************************************************
c This subroutine change the Jacobi coordinates and matrices of the tp
c  if a hierarchy change occurred within the bodies
c
c Inputs:
c    nbod          ==> Number of massive bodies
c    mass          ==> Masses of bodies
c    umatp, matp     ==> Conversion matrixes (Jacobi - Barycentric) for tp
c                       (old hierarchy)
c    oloctp ==> Link Bodies <--> tp's Orbits (old)
c    etatp       ==> Masses for tp's centers (old)
c    umat, mat     ==> Conversion matrixes (Jacobi - Barycentric) (new)
c    olocold, oloc ==> Link Bodies <--> Orbits (old and new)
c    eta, mu       ==> Masses for centers & sats (new)
c    xb, yb, zb    ==> Bodies position in Barycentric coord
c    xj, yj, zj    ==> Bodies position in Gen. Jacobi coord (old)
c    vxj, vyj, vzj ==> Bodies velocity in Gen. Jacobi coord (old)
c    axb, ayb, azb ==> Bodies acceleration in barycentric coord (new)
c    xjt, yjt, zjt ==> Tp position in Gen. Jacobi coord (old)
c    vxjt, vyjt, vzjt ==> Tp velocity in Gen. Jacobi coord (old)
c    dt            ==> Time step
c
c Outputs:
c    xj, yj, zj    ==> Bodies position in Gen. Jacobi coord (new)
c    vxj, vyj, vzj ==> Bodies velocity in Gen. Jacobi coord (new)
c    xjt, yjt, zjt ==> Tp position in Gen. Jacobi coord (new)
c    vxjt, vyjt, vzjt ==> Tp velocity in Gen. Jacobi coord (new)
c    axjt, ayjt, azjt ==> Tp acceleration in Gen. Jacobi coord (new)
c    umatp, matp     ==> Conversion matrixes (Jacobi - Barycentric) for tp
c                       (new)
c    oloctp      ==> Link Bodies <--> tp's Orbits (new)
c    etatp       ==> Masses for tp's centers (new)

      subroutine ce_change_hierarch_tp(i1st, nbod, mass, eta, mu,
     &     oloc, olocold, mat, umat, oloct, matp, umatp, etatp,
     &     xj, yj, zj, vxj, vyj, vzj, axb, ayb, azb, ir3j,
     &     xjt, yjt, zjt, vxjt, vyjt, vzjt, axjt, ayjt, azjt)

      include '../swift.inc'

c...  Inputs Only:
      integer nbod, oloc(NPLMAX, NPLMAX), olocold(nbod, nbod), i1st
      real*8 mass(nbod), eta(nbod), mu(nbod)
      real*8 mat(NPLMAX, NPLMAX), umat(NPLMAX,NPLMAX)
      real*8 axb(nbod), ayb(nbod), azb(nbod), ir3j(nbod)

c...  Inputs and Outputs:
      integer oloct(NPLMAX)
      real*8 etatp, umatp(NPLMAX), matp(NPLMAX)
      real*8 xjt, yjt, zjt, vxjt, vyjt, vzjt
      real*8 xj(nbod), yj(nbod), zj(nbod)
      real*8 vxj(nbod), vyj(nbod), vzj(nbod)

c...  Outputs Only:
      real*8 axjt, ayjt, azjt

c...  Internals:
      integer i, j, orbct(nbod)
      real*8 vsat, vcen, mtot
      real*8 xbt, ybt, zbt, vxbt, vybt, vzbt, axbt, aybt, azbt
      real*8 xb(NPLMAX), yb(NPLMAX), zb(NPLMAX)
      real*8 xjold(NPLMAX), yjold(NPLMAX), zjold(NPLMAX)
      real*8 vxjold(NPLMAX), vyjold(NPLMAX), vzjold(NPLMAX)
      real*8 vxb(NPLMAX), vyb(NPLMAX), vzb(NPLMAX)
      real*8 umatold(NPLMAX, NPLMAX), etaold(nbod), muold(nbod)
      logical sat, cen, ok

      save xb, yb, zb, vxb, vyb, vzb,
     &    xjold, yjold, zjold, vxjold, vyjold, vzjold

c----
c...  Executable code

      if (i1st.eq.0) then

         print*,"tp reorganization following hierarchy change"

         xjold(1:nbod) = xj(1:nbod)
         yjold(1:nbod) = yj(1:nbod)
         zjold(1:nbod) = zj(1:nbod)

         vxjold(1:nbod) = vxj(1:nbod)
         vyjold(1:nbod) = vyj(1:nbod)
         vzjold(1:nbod) = vzj(1:nbod)

         etaold = 0
         muold = 0
         do j = 2,nbod
            do i = 1,nbod
               if (olocold(j,i).eq.1) muold(j) = muold(j)+mass(i)
               if (olocold(j,i).eq.-1) etaold(j) = etaold(j)+mass(i)
            end do
         end do

         umatold = 0.0d0
         do i = 1,nbod
            umatold(i,1) = 1.0d0
         end do
         do j = 2,nbod
            vsat = etaold(j)/(muold(j)+etaold(j))
            vcen = -muold(j)/(muold(j)+etaold(j))
            do i = 1,nbod
               if (olocold(j,i).eq.1) umatold(i,j) = vsat
               if (olocold(j,i).eq.-1) umatold(i,j) = vcen
            end do
         end do

         call coord_g2b(nbod, umatold, mass, xj, yj, zj, vxj, vyj, vzj,
     &        xb, yb, zb, vxb, vyb, vzb)

         call coord_b2g(nbod, mat, mass, xb, yb, zb, vxb, vyb, vzb,
     &        xj, yj, zj, vxj, vyj, vzj)

      end if

      call coord_g2b_tp(nbod, umatp, xjold, yjold, zjold,
     &   vxjold, vyjold, vzjold, xjt, yjt, zjt, vxjt, vyjt, vzjt,
     &   xbt, ybt, zbt, vxbt, vybt, vzbt)

      call ce_change_orbit_tp(nbod, mass, eta, mu, oloc, orbct,
     &     xb, yb, zb, xbt, ybt, zbt)

      oloct = 0
      do j = 2, nbod
         sat = .true.
         cen = .true.
         do i = 1,nbod
            if (orbct(i).eq.-1) then
               sat = sat.and.(oloc(j,i).eq.1)
               cen = cen.and.(oloc(j,i).eq.-1)
            end if
         end do
         if (sat) oloct(j) = 1
         if (cen) oloct(j) = -1
      end do

      etatp = 0.0d0
      matp = 0.0d0
      umatp = 0.0d0
      do i = 1,nbod
         if (orbct(i).eq.-1) etatp = etatp + mass(i)
      end do
      do i = 1,nbod
         if (orbct(i).eq.-1) matp(i) = -mass(i)/etatp
      end do
      do i = 1,nbod
         if (orbct(i).eq.-1) then
            do j=1,nbod
               umatp(j) = umatp(j) - matp(i)*umat(i,j)
            end do
         end if
      end do

      call coord_vb2vg_tp(nbod, matp, xb, yb, zb,
     &     xbt, ybt, zbt, xjt, yjt, zjt)

      call coord_vb2vg_tp(nbod, matp, vxb, vyb, vzb,
     &     vxbt, vybt, vzbt, vxjt, vyjt, vzjt)

      call getacch_tp_hjs(nbod, mass, eta, mu, xj, yj, zj, xb, yb, zb,
     &     ir3j, oloct, etatp, xbt, ybt, zbt, xjt, yjt, zjt,
     &     axbt, aybt, azbt)

      call coord_vb2vg_tp(nbod, matp, axb, ayb, azb,
     &     axbt, aybt, azbt, axjt, ayjt, azjt)


      return
      end  ! ce_change_hierarch_tp
c-----------------------------------------------------------------------
