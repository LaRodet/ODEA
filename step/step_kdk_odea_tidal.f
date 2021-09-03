c***********************************************************************
c                            STEP_KDK_ODEA_TIDAL.F
c***********************************************************************
c This subroutine takes a step in generalized Jacobi coordinates (HJS)
c both massive and test particles
c
c Input:
c    i1st ==> = 0 if first step; = 1 not
c    time ==> current time
c    nbod ==> number of massive bodies
c    ntp ==>  number of massive bodies
c    oloc ==>  Link bodies - orbits
c    oloct  ==>  Link orbits - tp's
c    mass ==>  Masses of bodies
c    mat, umat ==>  Conversion matrixes for bodies
c    matp, umatp ==>  Conversion vectors for tp's
c    eta, mu ==> Masses for center & satellites for bodies
c    etatp ==> Masses for centers for tp's
c    xj, yj, zj ==> Initial position in gen. Jacobi coord
c    vxj, vyj, vzj ==> Initial velocity in gen. Jacobi coord
c    xjt, yjt, zjt ==> Initial tp position in gen. Jacobi coord
c    vxjt, vyjt, vzjt ==> Initial tp velocity in gen. Jacobi coord
c    istat ==> status of the test particles
c                 istat(i,1) = 0 ==> active:  = 1 not
c                 istat(i,2) = -1 ==> Danby did not work
c    rstat ==> status of the test particles
c    dt            ==>  Time step
c    iflgchk       ==>  Options flag (bit 6 and bit 7)
c    rtidal        ==> Radius of the massive bodies
c    atidal        ==> Dimensionless moment of inertia
c    sx, sy, sz    ==> Initial components of the spin
c    j2tidal       ==> Quadrupole moment
c
c Output:
c    xj, yj, zj ==> Final position in gen. Jacobi coord
c    vxj, vyj, vzj   ==>  final velocity in gen. Jacobi coord
c    xjt, yjt, zjt    ==>  final position in gen. Jacobi coord
c    vxjt, vyjt, vzjt ==>  final position in gen. Jacobi coord
c    sx, sy, sz    ==> Final components of the spin

      subroutine step_kdk_odea_tidal(i1st, time, nbod, ntp, oloc, oloct,
     &   mat, umat, matp, umatp, mass, eta, mu, etatp, xj, yj, zj, vxj,
     &   vyj, vzj, xjt, yjt, zjt, vxjt, vyjt, vzjt, istat, rstat, dt,
     &   iflgchk, atidal, rtidal, stidalx, stidaly, stidalz, j2tidal)

      include '../swift.inc'

c...  Inputs Only:
      integer nbod, ntp, i1st, iflgchk
      real*8 mass(nbod), dt, time
      real*8 mu(nbod), eta(nbod), etatp(ntp)
      real*8 umat(NPLMAX, NPLMAX), mat(NPLMAX, NPLMAX)
      real*8 umatp(NPLMAX, NTPMAX), matp(NPLMAX, NTPMAX)
      integer oloc(NPLMAX, NPLMAX), oloct(NPLMAX, NTPMAX)
      real*8 atidal(nbod), rtidal(nbod), j2tidal(nbod)

c...  Inputs and Outputs:
      integer istat(NTPMAX, NSTAT)
      real*8 rstat(NTPMAX, NSTATR)
      real*8 xj(nbod), yj(nbod), zj(nbod)
      real*8 vxj(nbod), vyj(nbod), vzj(nbod)
      real*8 xjt(ntp), yjt(ntp), zjt(ntp)
      real*8 vxjt(ntp), vyjt(ntp), vzjt(ntp)
      real*8 stidalx(nbod), stidaly(nbod), stidalz(nbod)

c...  Internals
      integer i1sttp, i, olocold(nbod,nbod)
      real*8 xbbeg(NPLMAX), ybbeg(NPLMAX), zbbeg(NPLMAX)
      real*8 axbbeg(NPLMAX), aybbeg(NPLMAX), azbbeg(NPLMAX)
      real*8 xjbeg(NPLMAX), yjbeg(NPLMAX), zjbeg(NPLMAX)
      real*8 xjend(NPLMAX), yjend(NPLMAX), zjend(NPLMAX)
      real*8 vxjbeg(NPLMAX),vyjbeg(NPLMAX),vzjbeg(NPLMAX)
      real*8 ir3jbeg(NPLMAX), ir3jend(NPLMAX)
      real*8 vxjh(NPLMAX), vyjh(NPLMAX), vzjh(NPLMAX)
      real*8 xbend(NPLMAX), ybend(NPLMAX), zbend(NPLMAX)
      real*8 axbend(NPLMAX), aybend(NPLMAX), azbend(NPLMAX)
      logical change

      save xbbeg, ybbeg, zbbeg, axbbeg, aybbeg, azbbeg, ir3jbeg
c----
c...  Executable code

      i1sttp = i1st

c...  remember the current position & velocities of the massive bodies
      do i=1, nbod
         xjbeg(i) = xj(i)
         yjbeg(i) = yj(i)
         zjbeg(i) = zj(i)
         vxjbeg(i) = vxj(i)
         vyjbeg(i) = vyj(i)
         vzjbeg(i) = vzj(i)
      enddo

      change = .False.

c...  first do the planets
      call step_kdk_pl_odea_tidal(i1st, time, nbod, oloc, umat, mat,
     &   mass, eta, mu, xj, yj, zj, vxj, vyj, vzj, xbbeg, ybbeg, zbbeg,
     &   axbbeg, aybbeg, azbbeg, ir3jbeg, vxjh, vyjh, vzjh, xbend,
     &   ybend, zbend, axbend, aybend, azbend, ir3jend, dt, olocold,
     &   change, iflgchk, atidal, rtidal, stidalx, stidaly, stidalz,
     &   j2tidal)

      if(ntp.ne.0) then

c...     now remember these positions
         do i=1,nbod
            xjend(i) = xj(i)
            yjend(i) = yj(i)
            zjend(i) = zj(i)
         enddo

c...     next the test particles
         call step_kdk_tp_odea(i1sttp, nbod, ntp, matp, umatp, oloct,
     &   mass,eta, mu, etatp, xjbeg, yjbeg, zjbeg, vxjbeg, vyjbeg,
     &   vzjbeg, ir3jbeg, xbbeg, ybbeg, zbbeg, axbbeg, aybbeg, azbbeg,
     &   vxjh, vyjh, vzjh, xjend, yjend, zjend, ir3jend, xbend, ybend,
     &   zbend, axbend, aybend, azbend, xjt, yjt, zjt, vxjt, vyjt, vzjt,
     &   istat, dt, oloc, olocold, mat, umat, change, time)

      endif

c...  store things for next step
      do i=1,nbod
         xbbeg(i) = xbend(i)
         ybbeg(i) = ybend(i)
         zbbeg(i) = zbend(i)
         axbbeg(i) = axbend(i)
         aybbeg(i) = aybend(i)
         azbbeg(i) = azbend(i)
         ir3jbeg(i) = ir3jend(i)
      end do

      return
      end   ! step_kdk_odea_tidal
c-----------------------------------------------------------------------
