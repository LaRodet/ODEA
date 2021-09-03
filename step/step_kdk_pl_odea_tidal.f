c***********************************************************************
c                        STEP_KDK_PL_ODEA_TIDAL.F
c***********************************************************************
c This subroutine takes a step in generalized Jacobi coord (HJS)
c Does a KICK then a DRIFT then a KICK.
c Only does massive particles.
c
c Input:
c    i1st          ==>  = 0 if first step; = 1 not
c    nbod          ==>  Number of massive bodies
c    mass          ==>  Masses of bodies
c    umat, mat     ==>  Conversion matrixes (Jacobi - Barycentric)
c    oloc          ==>  Link Bodies <--> Orbits
c    eta, mu       ==>  Masses for centers & sats
c    xj, yj, zj    ==>  Initial position in Gen. Jacobi coord
c    vxj, vyj, vzj ==>  Initial velocity in Gen. Jacobi coord
c    dt            ==>  Time step
c    rtidal        ==> Radius of the massive bodies
c    atidal        ==> Dimensionless moment of inertia
c    j2tidal       ==> Quadrupole moment
c
c Input & Output
c    xbbeg, ybbeg, zbbeg    ==>  Initial position in bary coord
c    axbbeg, aybbeg, azbbeg ==>  Initial accel. in bary coord
c    xj, yj, zj             ==>  Final position in Gen. Jacobi coord
c    vxj, vyj, vzj          ==>  Final velocity in Gen. Jacobi coord
c    ir3jbeg, ir3jend       ==>  Inverse radii^3 beg & end
c    xbend, ybend, zbend    ==>  Final position in bary coord
c    axbend, aybend, azbend ==>  Final accel in bary coord
c    vxjh, vyjh, vzjh       ==>  Middle velocity in Gen. Jacobi coord
c    sx, sy, sz             ==>  Components of the spin

      subroutine step_kdk_pl_odea_tidal(i1st, time, nbod, oloc, umat,
     &   mat, mass, eta, mu, xj, yj, zj, vxj, vyj, vzj, xbbeg, ybbeg,
     &   zbbeg, axbbeg, aybbeg, azbbeg, ir3jbeg, vxjh, vyjh, vzjh,
     &   xbend, ybend, zbend, axbend, aybend, azbend, ir3jend, dt,
     &   olocold, change, iflgchk, atidal, rtidal, stidalx, stidaly,
     &   stidalz, j2tidal)

      include '../swift.inc'

c...  Inputs Only:
      integer nbod, i1st, iflgchk
      real*8 mass(nbod), dt, time
      real*8 atidal(nbod), rtidal(nbod), j2tidal(nbod)

c...  Inputs and Outputs:
      integer oloc(NPLMAX, NPLMAX)
      real*8 eta(nbod), mu(nbod)
      real*8 umat(NPLMAX, NPLMAX), mat(NPLMAX, NPLMAX)
      real*8 xj(nbod), yj(nbod), zj(nbod)
      real*8 vxj(nbod), vyj(nbod), vzj(nbod)
      real*8 xbbeg(nbod), ybbeg(nbod), zbbeg(nbod)
      real*8 stidalx(nbod), stidaly(nbod), stidalz(nbod)

c...  Outputs Only:
      integer olocold(nbod, nbod)
      logical change
      real*8 ir3jbeg(nbod), ir3jend(nbod)
      real*8 xbend(nbod), ybend(nbod), zbend(nbod)
      real*8 vxjh(nbod), vyjh(nbod), vzjh(nbod)
      real*8 axbbeg(nbod),aybbeg(nbod),azbbeg(nbod)
      real*8 axbend(nbod),aybend(nbod),azbend(nbod)

c...  Internals:
      integer i, j, k
      logical ok, checkchange
      real*8 dth
      real*8 axj(NPLMAX), ayj(NPLMAX), azj(NPLMAX)
      real*8 vxb(NPLMAX), vyb(NPLMAX), vzb(NPLMAX)

      Save axj, ayj, azj, checkchange

c----
c...  Executable code

      if (i1st.eq.0) then
c...     Compute barycentric coords
         call  coord_g2b(nbod, umat, mass, xj, yj, zj, vxj, vyj, vzj,
     &      xbbeg, ybbeg, zbbeg, vxb, vyb, vzb)
c...     Get the Jacobi accelerations if first time step
         call getaccj_hjs(nbod, oloc, mass, eta, mu, xj, yj, zj,
     &      xbbeg, ybbeg, zbbeg, ir3jbeg, axj, ayj, azj)
         call tidal_quadrupole(nbod, oloc, mass, eta, mu, xj, yj, zj,
     &      axj, ayj, azj, rtidal, stidalx, stidaly, stidalz, j2tidal)
         call coord_g2b(nbod, umat, mass, vxj, vyj, vzj, axj, ayj, azj,
     &      vxb, vyb, vzb, axbbeg, aybbeg, azbbeg)
         checkchange = btest(iflgchk,6)   ! Initial check
         i1st = 1               ! turn this off
      endif

      if (checkchange) then
c         print*,'t=',time,' considering change'

         call ce_change_hierarch(time, nbod, oloc, umat, mat, mass, eta,
     &   mu, xj, yj, zj, vxj, vyj, vzj, axj(1:nbod), ayj(1:nbod),
     &   azj(1:nbod), xbbeg, ybbeg, zbbeg, axbbeg, aybbeg, azbbeg,
     &   ir3jbeg, dt, olocold, change, iflgchk)

      end if

      dth = 0.5d0*dt

c...  Apply a Jacobi kick for a half dt
      call kickvh(nbod, vxj, vyj, vzj, axj, ayj, azj, dth)
      call tidal_quadrupole_spin_evolve(nbod, oloc, mass, eta, mu,
     &     xj, yj, zj, atidal, rtidal, stidalx, stidaly, stidalz,
     &     j2tidal, dth)

c...  Drift in Jacobi coords for the full step
      call drift_hjs(nbod, mass, eta, mu, xj, yj, zj, vxj, vyj, vzj, dt)

c...  After drift, compute bary. xb and vb for acceleration calculations
      call coord_g2b(nbod, umat, mass, xj, yj, zj, vxj, vyj, vzj,
     &     xbend, ybend, zbend, vxb, vyb, vzb)

c...  Save Jacobi velocities at middle point (for tp's after)
      do i=1, nbod
         vxjh(i) = vxj(i)
         vyjh(i) = vyj(i)
         vzjh(i) = vzj(i)
      end do

      call getaccj_hjs(nbod, oloc, mass, eta, mu, xj, yj, zj,
     &     xbend, ybend, zbend, ir3jend, axj, ayj, azj)
      call tidal_quadrupole(nbod, oloc, mass, eta, mu, xj, yj, zj,
     &      axj, ayj, azj, rtidal, stidalx, stidaly, stidalz, j2tidal)
      call ce_crit(nbod, xj, yj, zj, mass, mu, eta, oloc,
     &     axj, ayj, azj, checkchange)

      call coord_g2b(nbod, umat, mass, vxj, vyj, vzj, axj, ayj, azj,
     &     vxb, vyb, vzb, axbend, aybend, azbend)

c...  Apply a Jacobi kick for a half dt
      call kickvh(nbod, vxj, vyj, vzj, axj, ayj, azj, dth)
      call tidal_quadrupole_spin_evolve(nbod, oloc, mass, eta, mu,
     &     xj, yj, zj, atidal, rtidal, stidalx, stidaly, stidalz,
     &     j2tidal, dth)

      return
      end   ! step_kdk_pl_odea
c-----------------------------------------------------------------------
