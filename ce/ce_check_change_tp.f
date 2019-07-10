c***********************************************************************
c                            CE_CHECK_CHANGE_TP.F
c***********************************************************************
c This subroutine check is a tp should change orbit
c
c Inputs:
c    nbod          ==> Number of massive bodies
c    mass          ==> Masses of bodies
c    umatp, matp   ==> Conversion matrixes (Jacobi - Barycentric) for tp
c    oloctp        ==> Link Bodies <--> tp's Orbits
c    etatp         ==> Masses for tp's centers
c    umat, mat     ==> Conversion matrixes (Jacobi - Barycentric)
c    oloc          ==> Link Bodies <--> Orbits
c    eta, mu       ==> Masses for centers & sats
c    xj, yj, zj    ==> Bodies position in Gen. Jacobi coord
c    vxj, vyj, vzj ==> Bodies velocity in Gen. Jacobi coord
c    axb, ayb, azb ==> Bodies acceleration in barycentric coord
c    xjt, yjt, zjt ==> Tp position in Gen. Jacobi coord
c    vxjt, vyjt, vzjt ==> Tp velocity in Gen. Jacobi coord
c    time            ==> Time
c
c Outputs:
c    xjt, yjt, zjt ==> Tp position in Gen. Jacobi coord (new)
c    vxjt, vyjt, vzjt ==> Tp velocity in Gen. Jacobi coord (new)
c    axjt, ayjt, azjt ==> Tp acceleration in Gen. Jacobi coord (new)
c    umatp, matp     ==> Conversion matrixes (Jacobi - Barycentric) for tp
c                       (new)
c    oloctp      ==> Link Bodies <--> tp's Orbits (new)
c    etatp       ==> Masses for tp's centers (new)

      subroutine ce_check_change_tp(num, nbod, mass, eta, mu, oloc,
     &   umat, oloct, matp, umatp, etatp, xj, yj, zj, vxj, vyj, vzj,
     &   axb, ayb, azb, ir3j, xjt, yjt, zjt, vxjt, vyjt, vzjt,
     &   axjt, ayjt, azjt, time)

      include '../swift.inc'

c...  Inputs Only:
      integer nbod, oloc(NPLMAX, NPLMAX), num
      real*8 mass(nbod), eta(nbod), mu(nbod), time
      real*8 umat(NPLMAX, NPLMAX)
      real*8 xj(nbod), yj(nbod), zj(nbod)
      real*8 vxj(nbod), vyj(nbod), vzj(nbod)
      real*8 axb(nbod), ayb(nbod), azb(nbod), ir3j(nbod)

c...  Inputs and Outputs
      integer oloct(NPLMAX)
      real*8 umatp(NPLMAX), matp(NPLMAX), etatp
      real*8 xjt, yjt, zjt, vxjt, vyjt, vzjt, axjt, ayjt, azjt

c...  Internals
      integer i, j
      real*8 xbt, ybt, zbt, vxbt, vybt, vzbt, axbt, aybt, azbt
      real*8 xb(nbod), yb(nbod), zb(nbod)
      real*8 vxb(nbod), vyb(nbod), vzb(nbod)
      integer orbctold(nbod), orbct(nbod)
      logical sat, cen

c----
c...  Executable code

      call coord_g2b(nbod, umat, mass, xj, yj, zj, vxj, vyj, vzj,
     &     xb, yb, zb, vxb, vyb, vzb)

      do i=1,nbod
         if (matp(i).ne.0.d0) then
            orbctold(i)=-1
         else
            orbctold(i)=0
         endif
      end do

      call coord_g2b_tp(nbod, umatp, xj, yj, zj, vxj, vyj, vzj, xjt,
     &   yjt, zjt, vxjt, vyjt, vzjt, xbt, ybt, zbt, vxbt, vybt, vzbt)

      call ce_change_orbit_tp(nbod, mass, eta, mu, oloc, orbct,
     &   xb, yb, zb, xbt, ybt, zbt)

      if (any(orbctold.ne.orbct)) then

c         print*, "t", time, "tp", num, "changes orbit"
c         print*, "Centers:", orbct(1:nbod)

         oloct = 0
         do j = 2,nbod
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
     &      xbt, ybt, zbt, xjt, yjt, zjt)

         call coord_vb2vg_tp(nbod,matp,vxb,vyb,vzb,
     &      vxbt,vybt,vzbt,vxjt,vyjt,vzjt)

         call getacch_tp_hjs(nbod, mass, eta, mu, xj, yj, zj, xb, yb,
     &      zb, ir3j, oloct, etatp,
     &      xbt, ybt, zbt, xjt, yjt, zjt, axbt, aybt, azbt)

         call coord_vb2vg_tp(nbod, matp, axb, ayb, azb,
     &      axbt, aybt, azbt, axjt, ayjt, azjt)

      end if

      return

      end ! ce_check_change_tp
c-----------------------------------------------------------------------
