c***********************************************************************
c                           TIDAL_DISSIPATION.F
c***********************************************************************
c Given the cartesian position and velocity of an orbit and the tidal
c parameter, compute the secular tidal correction.
c
c ALGORITHM: Correia et al. 2011

      subroutine tidal_dissipation(nbod, oloc, mass, eta, mu,
     &     xj, yj, zj, vxj, vyj, vzj, atidal, qtidal, rtidal, sx, sy,
     &     sz, dt)

      include '../swift.inc'

c...  Inputs Only:
      real*8 mass(NPLMAX), atidal(NPLMAX), qtidal(NPLMAX)
      real*8 rtidal(NPLMAX)
      real*8 mu(NPLMAX), eta(NPLMAX), dt, mm, mc, mp
      integer oloc(NPLMAX,NPLMAX)
      integer nbod

c...  Inputs and Outputs
      real*8 xj(NPLMAX), yj(NPLMAX), zj(NPLMAX)
      real*8 vxj(NPLMAX), vyj(NPLMAX), vzj(NPLMAX)
      real*8 sx(NPLMAX), sy(NPLMAX), sz(NPLMAX), s

c...  Internals:
      real*8 a, e, inc, capom, omega, capm, gm, da, de, beta
      real*8 dsx, dsy, dsz, dsfact, defact
      real*8 f1, f2, f3, f4, f5, e2, e4, e6, e8, r
      real*8 kx, ky, kz, k, ct, ex, ey, ez, es
      real*8 dex, dey, dez, dkx, dky, dkz
      integer i, j, n, ialpha, itidal, icomp
      logical tidal

c----
c...  Executable code
      do j=2,nbod
        n = 0
        tidal = .false.
        do i=1,nbod
          if (oloc(j,i).ne.0) then
            n = n+1
            if (qtidal(i).gt.0d0) then
              tidal = .true.
              itidal = i
            else
              icomp = i
            end if
          end if
        end do
        if ((n.eq.2) .and. tidal) then
          gm = eta(j)+mu(j)
          call orbel_xv2el(xj(j), yj(j), zj(j), vxj(j), vyj(j), vzj(j),
     &         gm, ialpha, a, e, inc, capom, omega, capm)
          if (ialpha.eq.-1) then
            beta = eta(j)/(4*PI*PI)*mu(j)/(eta(j)+mu(j))
            k = sqrt(gm*a*(1-e*e))*beta
            kx = k*sin(inc)*sin(capom)
            ky = -k*sin(inc)*cos(capom)
            kz = k*cos(inc)

            s = sqrt(sx(itidal)*sx(itidal) + sy(itidal)*sy(itidal) +
     &          sz(itidal)*sz(itidal))
            ct = (kx*sx(itidal) + ky*sy(itidal) + kz*sz(itidal))/(s*k)

            e2 = e*e
            e4 = e2*e2
            e6 = e2*e4
            e8 = e4*e4
            ex =e*(cos(capom)*cos(omega)-cos(inc)*sin(capom)*sin(omega))
            ey =e*(sin(capom)*cos(omega)+cos(inc)*cos(capom)*sin(omega))
            ez = e*sin(inc)*sin(omega)
            es = (ex*sx(itidal)+ey*sy(itidal)+ez*sz(itidal))/s
            f1 = (1d0+3d0*e2+3d0*e4/8d0)/((1d0-e2)**4.5d0)
            f2 = (1d0+7.5d0*e2+45d0*e4/8d0+5d0*e6/16d0)/((1d0-e2)**6d0)
            f3 = (1d0+15.5d0*e2+255d0*e4/8d0+185d0*e6/16d0+25d0*e8/64d0)
     &           /((1d0-e2)**7.5d0)
            f4 = (1d0+1.5d0*e2+e4/8d0)/((1d0-e2)**5d0)
            f5 = (1+15d0*e2/4d0+15d0*e4/8d0+5*e6/64d0)/((1-e2)**6.5d0)

            mm = sqrt(gm/(a**3))
            r = rtidal(itidal)
            mc = mass(icomp)
            mp = mass(itidal)

            dsfact = 9d0/2d0*dt*mm/(qtidal(itidal)*atidal(itidal))*
     &         ((r/a)**3)*(mc**2)/(mp*(mp+mc))
            dsx =  dsfact*(f4*sqrt(1-e2)*0.5d0*s*(sx(itidal)/s-ct*kx/k)-
     &             f1*sx(itidal)+f2*mm*kx/k+
     &             es*(6d0+e2)*s*ex/((4d0*(1-e2)**4.5d0)))
            dsy =  dsfact*(f4*sqrt(1-e2)*0.5d0*s*(sy(itidal)/s-ct*ky/k)-
     &             f1*sy(itidal)+f2*mm*ky/k
     &             +es*(6d0+e2)*s*ey/((4d0*(1-e2)**4.5d0)))
            dsz =  dsfact*(f4*sqrt(1-e2)*0.5d0*s*(sz(itidal)/s-ct*kz/k)-
     &             f1*sz(itidal)+f2*mm*kz/k+
     &             es*(6d0+e2)*s*ez/((4d0*(1-e2)**4.5d0)))

            dkx = -dsx*(mp/(4*PI*PI))*(r**2)*atidal(itidal)
            dky = -dsy*(mp/(4*PI*PI))*(r**2)*atidal(itidal)
            dkz = -dsz*(mp/(4*PI*PI))*(r**2)*atidal(itidal)

            kx = kx + dkx
            ky = ky + dky
            kz = kz + dkz
            k = sqrt(kx*kx + ky*ky + kz*kz)

            defact = -9d0/2d0*dt/qtidal(itidal)*((r/a)**5)*(mc/mp)
            dex = defact*(f4*0.5d0*s*es*kx/k-(5.5d0*f4*ct*s-9*f5*mm)*ex)
            dey = defact*(f4*0.5d0*s*es*ky/k-(5.5d0*f4*ct*s-9*f5*mm)*ey)
            dez = defact*(f4*0.5d0*s*es*kz/k-(5.5d0*f4*ct*s-9*f5*mm)*ez)

            ex = ex+dex
            ey = ey+dey
            ez = ez+dez
            e = sqrt(ex*ex+ey*ey+ez*ez)

            a = (k/beta)**2 / (gm*(1-e*e))

            inc = acos(kz/k)
            if (inc.lt.TINY) then
              omega = mod(atan2(ey, ex),(2*PI))
            else
              omega = mod(atan2(ez, (kx*ey-ky*ex)/k),(2*PI))
            end if
            capom = mod(atan2(kx, -ky),(2*PI))

c            da = dt*9d0*a/qtidal(itidal)* ((r/a)**5)*(mc/mp)*
c     &         (f2*ct*s-mm*f3)
c            de = 81d0/2d0*dt/qtidal(itidal)*((r/a)**5)*(mc/mp)*
c     &         (11d0*f4*ct*s/18d0-mm*f5)*e

            sx(itidal) = sx(itidal)+dsx
            sy(itidal) = sy(itidal)+dsy
            sz(itidal) = sz(itidal)+dsz

c            write(*,*)a, e, inc, capom, omega

            call orbel_el2xv(gm, ialpha, a, e, inc, capom, omega, capm,
     &          xj(j), yj(j), zj(j), vxj(j), vyj(j), vzj(j))
          end if
        end if
      end do

      return
      end    ! tidal_dissipation
c-----------------------------------------------------------------------
