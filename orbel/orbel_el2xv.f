c***********************************************************************
c                          ORBEL_EL2XV.F
c***********************************************************************
c PURPOSE: To compute cartesian positions and velocities given
c  central mass, ialpha ( = +1 for hyp., 0 for para. and
c  -1 for ellipse), and orbital elements.
c
c Inputs:
c    gm       ==> G times central mass
c    ialpha   ==> conic section type
c    a        ==> semi-major axis or pericentric distance if a parabola
c    e        ==> eccentricity
c    inc      ==> inclination
c    capom    ==> longitude of ascending node
c    omega    ==> argument of perihelion
c    capm     ==> mean anomoly
c
c Outputs:
c    x, y, z    ==>  position of object
c    vx, vy, vz ==>  velocity of object
c
c ALGORITHM:  See Fitzpatrick "Principles of Cel. Mech."
c REMARKS: All angles are in RADIANS

      subroutine orbel_el2xv(gm, ialpha, a, e, inc, capom, omega, capm,
     &   x, y, z, vx, vy, vz)

      include '../swift.inc'

c...  Inputs Only:
      integer ialpha
      real*8 gm, a, e, inc, capom, omega, capm

c...  Outputs:
      real*8 x, y, z, vx, vy, vz

c...  Internals:
      real*8 cape, capf, zpara, em1
      real*8 sp, cp, so, co, si, ci
      real*8 d11, d12, d13, d21, d22, d23
      real*8 scap, ccap, shcap, chcap
      real*8 sqe, sqgma, xfac1, xfac2, ri, vfac1, vfac2
      real*8 orbel_ehybrid, orbel_fhybrid, orbel_zget

c----
c...  Executable code

      if(e.lt.0.0) then
         write(*,*) ' ERROR in orbel_el2xv: e<0, setting e=0!!1'
         e = 0.0
      endif

c...  Check for inconsistencies between ialpha and e
      em1 = e - 1.d0
      if(
     &     ((ialpha.eq.0) .and. (abs(em1).gt.TINY))  .or.
     &     ((ialpha.lt.0) .and. (e.gt.1.0d0))  .or.
     &     ((ialpha.gt.0) .and. (e.lt.1.0d0)) )  then
         write(*,*) 'ERROR in orbel_el2xv: ialpha and e inconsistent'
         write(*,*) '                       ialpha = ',ialpha
         write(*,*) '                            e = ',e
      endif

c...  Generate rotation matrices (on p. 42 of Fitzpatrick)
      call orbel_scget(omega, sp, cp)
      call orbel_scget(capom, so, co)
      call orbel_scget(inc, si, ci)
      d11 = cp*co - sp*so*ci
      d12 = cp*so + sp*co*ci
      d13 = sp*si
      d21 = -sp*co - cp*so*ci
      d22 = -sp*so + cp*co*ci
      d23 = cp*si

c...  Get the other quantities depending on orbit type ( i.e. IALPHA)
      if (ialpha .eq. -1) then
         cape = orbel_ehybrid(e,capm)
         call orbel_scget(cape,scap,ccap)
         sqe = sqrt(1.d0 -e*e)
         sqgma = sqrt(gm*a)
         xfac1 = a*(ccap - e)
         xfac2 = a*sqe*scap
         ri = 1.d0/(a*(1.d0 - e*ccap))
         vfac1 = -ri * sqgma * scap
         vfac2 = ri * sqgma * sqe * ccap
      endif

      if (ialpha .eq. +1) then
         capf = orbel_fhybrid(e,capm)
         call orbel_schget(capf, shcap, chcap)
         sqe = sqrt(e*e - 1.d0 )
         sqgma = sqrt(gm*a)
         xfac1 = a*(e - chcap)
         xfac2 = a*sqe*shcap
         ri = 1.d0/(a*(e*chcap - 1.d0))
         vfac1 = -ri * sqgma * shcap
         vfac2 = ri * sqgma * sqe * chcap
      endif

      if (ialpha .eq. 0) then
         zpara = orbel_zget(capm)
         sqgma = sqrt(2.d0*gm*a)
         xfac1 = a*(1.d0 - zpara*zpara)
         xfac2 = 2.d0*a*zpara
         ri = 1.d0/(a*(1.d0 + zpara*zpara))
         vfac1 = -ri * sqgma * zpara
         vfac2 = ri * sqgma
      endif

      x =  d11*xfac1 + d21*xfac2
      y =  d12*xfac1 + d22*xfac2
      z =  d13*xfac1 + d23*xfac2
      vx = d11*vfac1 + d21*vfac2
      vy = d12*vfac1 + d22*vfac2
      vz = d13*vfac1 + d23*vfac2

      return
      end    ! orbel_el2xv
c-----------------------------------------------------------------------
