c***********************************************************************
c                           TIDAL_DISSIPATION.F
c***********************************************************************
c Given the cartesian position and velocity of an orbit,
c       compute the osculating orbital elements.
c
c Inputs:
c    x, y, z    ==> position of object
c    vx, vy, vz ==> velocity of object
c    gmsum      ==> G * total mass
c
c Outputs:
c    ialpha     ==> conic section type
C    a          ==> semi-major axis or pericentric distance if a parabola
c    e          ==> eccentricity
C    inc        ==> inclination
C    capom      ==> longitude of ascending node
C    omega      ==> argument of perihelion
C    capm       ==> mean anomaly
c
c ALGORITHM: See e.g. p.70 of Fitzpatrick's "Priciples of Cel. Mech."
c REMARKS:  If the inclination INC is less than TINY, we arbitrarily
c            choose the longitude of the ascending node LGNODE to be 0.0
c            (so the ascending node is then along the X axis).
c           If the  eccentricity E is less than SQRT(TINY),
c            we arbitrarily choose the argument of perihelion to be 0.

      subroutine tidal_dissipation(nbod, oloc, mass, eta, mu,
     &     xj, yj, zj, vxj, vyj, vzj, atidal, qtidal, rtidal, stidal,
     &     dt)

      include '../swift.inc'

c...  Inputs Only:
      real*8 mass(NPLMAX), atidal(NPLMAX), qtidal(NPLMAX)
      real*8 rtidal(NPLMAX)
      real*8 mu(NPLMAX), eta(NPLMAX), dt, meanmotion
      integer oloc(NPLMAX,NPLMAX)
      integer nbod

c...  Inputs and Outputs
      real*8 xj(NPLMAX), yj(NPLMAX), zj(NPLMAX)
      real*8 vxj(NPLMAX), vyj(NPLMAX), vzj(NPLMAX)
      real*8 stidal(NPLMAX)

c...  Internals:
      real*8 a, e, inc, capom, omega, capm, gm, da, ds
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
            meanmotion = sqrt(gm/(a**3))
            da = dt*9d0*a*meanmotion/qtidal(itidal)*
     &         ((rtidal(itidal)/a)**5)*mass(icomp)/mass(itidal)*
     &         (stidal(itidal)-meanmotion)/meanmotion
            ds = -9d0/2d0*dt*(meanmotion**2)/(qtidal(itidal)*
     &         atidal(itidal))*((rtidal(itidal)/a)**3)*mass(icomp)
     &         /mass(itidal)*(stidal(itidal)-meanmotion)/meanmotion
            a = a+da
            stidal(itidal) = stidal(itidal)+ds
            call orbel_el2xv(gm, ialpha, a, e, inc, capom, omega, capm,
     &          xj(j), yj(j), zj(j), vxj(j), vyj(j), vzj(j))
          end if
        end if
      end do

      return
      end    ! tidal_dissipation
c-----------------------------------------------------------------------
