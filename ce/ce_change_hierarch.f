c***********************************************************************
c                            CE_CHANGE_HIERARCH.F
c***********************************************************************
c This subroutine finds the best hierarchy according to an acceleration
c criterion and update the coordinates and matrices if needed.
c Only does massive particles
c
c Inputs:
c    nbod          ==> Number of massive bodies
c    mass          ==> Masses of bodies
c    umat, mat     ==> Conversion matrixes (Jacobi - Barycentric)
c    oloc          ==> Link Bodies <--> Orbits
c    eta, mu       ==> Masses for centers & sats
c    xb, yb, zb    ==> Position in Barycentric coord
c    xj, yj, zj    ==> Position in Gen. Jacobi coord (current hierarchy)
c    vxj, vyj, vzj ==> Velocity in Gen. Jacobi coord (current hierarchy)
c    dt            ==> Time step
c
c Outputs:
c    xj, yj, zj    ==> Position in Gen. Jacobi coord (new hierarchy)
c    vxj, vyj, vzj ==> Velocity in Gen. Jacobi coord (new hierarchy)
c    ir3jbeg, ir3jend ==>  Inverse radii^3 beg & end
c    axb, ayb, azb ==> Acceleration in barycentric coord (new hierarchy)
c    olocold, oloc ==>  Link Bodies <--> Orbits, old and new
c    change ==> True if a hierarchy change occured


      subroutine ce_change_hierarch(time, nbod, oloc, umat, mat, mass,
     &   eta, mu, xj, yj, zj, vxj, vyj, vzj, axj, ayj, azj, xb, yb, zb,
     &   axb, ayb, azb, ir3j, dt, olocold, change, iflgchk)

      include '../swift.inc'

c...  Inputs Only:
      integer nbod, iflgchk
      real*8 xb(nbod), yb(nbod), zb(nbod), time

c...  Inputs and Outputs:
      integer oloc(NPLMAX,NPLMAX)
      real*8 mass(nbod), eta(nbod), mu(nbod)
      real*8 umat(NPLMAX,NPLMAX), mat(NPLMAX,NPLMAX)
      real*8 xj(nbod), yj(nbod), zj(nbod)
      real*8 vxj(nbod), vyj(nbod), vzj(nbod)

c...  Outputs Only
      real*8 dt
      integer olocold(nbod,nbod)
      logical change
      real*8 axb(nbod), ayb(nbod), azb(nbod)
      real*8 axj(nbod), ayj(nbod), azj(nbod)
      real*8 ir3j(nbod)

c...  Internals:
      integer i, j, k, l, imax, jmax
      integer olocnew(nbod,nbod)
      real*8 vxb(nbod), vyb(nbod), vzb(nbod)
      real*8 a(nbod,nbod), maxa
      real*8 rij2, mtot, etanew(nbod), munew(nbod)
      real*8 xg(nbod), yg(nbod), zg(nbod)

c----
c...  Executable code

      a = reshape((/ (-1, i=1,nbod*nbod) /), shape(a))

      do i=1,nbod-1
         do j=i+1,nbod
            rij2 = (xb(i)-xb(j))**2+(yb(i)-yb(j))**2+(zb(i)-zb(j))**2
            a(i,j) = (mass(j)+mass(i))/rij2
         end do
      end do

      xg = 0.0d0
      yg = 0.0d0
      zg = 0.0d0
      xg(1) = xj(1)
      yg(1) = yj(1)
      zg(1) = zj(1)

      olocold(1:nbod,1:nbod) = oloc(1:nbod,1:nbod)

      olocnew = 0
      munew = 0.0d0
      etanew = 0.0d0
      do k=2,nbod

         maxa = -1.
         imax = -1
         jmax = -1
         do i=1,nbod-1
            do j=i+1,nbod
               if (a(i,j).gt.maxa) then
                  maxa = a(i,j)
                  imax = i
                  jmax = j
               end if
            end do
         end do

         olocnew(k,imax) = -1
         etanew(k) = mass(imax)
         olocnew(k,jmax) = 1
         munew(k) = mass(jmax)

         do l=2,k-1
            if (olocnew(l,imax).ne.0) then
               do i=1,nbod
                  if ((olocnew(l,i).ne.0).and.(olocnew(k,i).eq.0)) then
                     olocnew(k,i) = -1
                     etanew(k) = etanew(k) + mass(i)
                  end if
               end do
            end if
            if (olocnew(l,jmax).ne.0) then
               do j=1,nbod
                  if ((olocnew(l,j).ne.0).and.(olocnew(k,j).eq.0)) then
                     olocnew(k,j) = 1
                     munew(k) = munew(k) + mass(j)
                  end if
               end do
            end if
         end do

         do i=1,nbod
            if (olocnew(k,i).ne.0) then
               xg(k) = xg(k) + mass(i)*xb(i)/(munew(k)+etanew(k))
               yg(k) = yg(k) + mass(i)*yb(i)/(munew(k)+etanew(k))
               zg(k) = zg(k) + mass(i)*zb(i)/(munew(k)+etanew(k))
            end if
         end do

         do i=1,nbod-1
            do j=i+1,nbod
               if ((olocnew(k,i).ne.0) .and. (olocnew(k,j).ne.0)) then
                  a(i,j) = -1
               elseif ((olocnew(k,i).eq.0) .and. (olocnew(k,j).ne.0))
     &                 then
                  rij2 = (xb(i)-xg(k))**2+(yb(i)-yg(k))**2+
     &                 (zb(i)-zg(k))**2
                  a(i,j) = (munew(k)+etanew(k)+mass(i))/rij2
               elseif ((olocnew(k,i).ne.0) .and. (olocnew(k,j).eq.0))
     &                 then
                  rij2 = (xg(k)-xb(j))**2+(yg(k)-yb(j))**2+
     &                 (zg(k)-zb(j))**2
                  a(i,j) =(munew(k)+etanew(k)+mass(j))/rij2
               end if
            end do
         end do

      end do

      change = .False.
      k = 1
      do while ((.not.change).and.(k.le.nbod))
         l = 1
         do while(any(oloc(k,1:nbod).ne.olocnew(l,1:nbod))
     &        .and.any(-oloc(k,1:nbod).ne.olocnew(l,1:nbod))
     &        .and.(l.le.nbod))
            l = l+1
         end do
         if (l.eq.nbod+1) change = .True.
         k = k+1
      end do

      if (change) then

         call coord_g2b(nbod, umat, mass, vxj, vyj, vzj, axj, ayj, azj,
     &      vxb, vyb, vzb, axb, ayb, azb)

         oloc(1:nbod,1:nbod) = olocnew(1:nbod,1:nbod)
         eta(1:nbod) = etanew(1:nbod)
         mu(1:nbod) = munew(1:nbod)

         mtot = 0.0d0
         do i = 1,nbod
            mtot = mtot+mass(i)
         end do
         mat = 0.0d0
         do i = 1,nbod
            mat(1,i) = mass(i)/mtot
         end do
         do j = 2,nbod
            do i = 1,nbod
               if (oloc(j,i).eq.1) mat(j,i) = mass(i)/mu(j)
               if (oloc(j,i).eq.-1) mat(j,i) = -mass(i)/eta(j)
            end do
         end do

         umat = 0.0d0
         do i = 1,nbod
            umat(i,1) = 1.0d0
         end do
         do j = 2,nbod
            do i = 1,nbod
               if (oloc(j,i).eq.1) umat(i,j) = eta(j)/(mu(j)+eta(j))
               if (oloc(j,i).eq.-1) umat(i,j) = -mu(j)/(mu(j)+eta(j))
            end do
         end do

         call coord_b2g(nbod, mat, mass, xb, yb, zb, vxb, vyb, vzb,
     &      xj, yj, zj, vxj, vyj, vzj)

         call getaccj_hjs(nbod, oloc, mass, eta, mu, xj, yj, zj,
     &      xb, yb, zb, ir3j, axj, ayj, azj)

         call coord_g2b(nbod, umat, mass, vxj, vyj, vzj, axj, ayj, azj,
     &      vxb, vyb, vzb, axb, ayb, azb)

         if (btest(iflgchk,7))
     &      call ce_change_dt(nbod, eta, mu, xj, yj, zj, vxj, vyj,
     &         vzj, dt)

         print*, 'Change of hierarchy at t=', time
         print*, 'New time step is', dt

      end if

      return
      end  ! ce_change_hierarch
c-----------------------------------------------------------------------
