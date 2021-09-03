c***********************************************************************
c                           TIDAL_QUADRUPOLE.F
c***********************************************************************
c Given the cartesian position and velocity of an orbit and the tidal
c parameter, compute the force associated with the quadrupole.
c
c ALGORITHM: Hong et al. 2018

      subroutine tidal_quadrupole(nbod, oloc, mass, eta, mu,
     &     xj, yj, zj, axj, ayj, azj, rtidal, sx, sy, sz, j2tidal)

      include '../swift.inc'

c...  Inputs Only:
      real*8 mass(NPLMAX)
      real*8 rtidal(NPLMAX), j2tidal(NPLMAX)
      real*8 mu(NPLMAX), eta(NPLMAX), mm, mc, mp
      integer oloc(NPLMAX,NPLMAX)
      integer nbod
      real*8 xj(NPLMAX), yj(NPLMAX), zj(NPLMAX)
      real*8 sx(NPLMAX), sy(NPLMAX), sz(NPLMAX), s

c...  Inputs and Outputs
      real*8 axj(NPLMAX), ayj(NPLMAX), azj(NPLMAX)

c...  Internals:
      real*8 gm, r2, r, j2, rs
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
            if (j2tidal(i).gt.0d0) then
              tidal = .true.
              itidal = i
            else
              icomp = i
            end if
          end if
        end do
        if ((n.eq.2) .and. tidal) then
          gm = eta(j)+mu(j)
          r2 = xj(j)*xj(j) + yj(j)*yj(j) + zj(j)*zj(j)
          r = sqrt(r2)
          s = sqrt(sx(itidal)*sx(itidal) + sy(itidal)*sy(itidal) +
     &        sz(itidal)*sz(itidal))
          rs = xj(j)*sx(itidal)+yj(j)*sy(itidal)+zj(j)*sz(itidal)
          rs = rs/(s*r)
          j2 = j2tidal(itidal)

          axj(j) = axj(j) + gm/r2*1.5d0*j2*(rtidal(itidal)/r)**2*
     &        ((5d0*rs*rs-1d0)*xj(j)/r - 2d0*rs*sx(itidal)/s)
          ayj(j) = ayj(j) + gm/r2*1.5d0*j2*(rtidal(itidal)/r)**2*
     &        ((5d0*rs*rs-1d0)*yj(j)/r - 2d0*rs*sy(itidal)/s)
          azj(j) = azj(j) + gm/r2*1.5d0*j2*(rtidal(itidal)/r)**2*
     &        ((5d0*rs*rs-1d0)*zj(j)/r - 2d0*rs*sz(itidal)/s)
        end if
      end do

      return
      end    ! tidal_quadrupole
c-----------------------------------------------------------------------
