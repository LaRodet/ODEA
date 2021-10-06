c***********************************************************************
c                TIDAL_QUADRUPOLE_SPIN_EVOLVE.F
c***********************************************************************
c Given the cartesian position and velocity of an orbit and the tidal
c parameter, compute the force associated with the quadrupole.
c
c ALGORITHM: Hong et al. 2018

      subroutine tidal_quadrupole_spin_evolve(nbod, oloc, mass, eta, mu,
     &     xj, yj, zj, atidal, rtidal, sx, sy, sz, j2tidal, dt)

      include '../swift.inc'

c...  Inputs Only:
      real*8 mass(NPLMAX), atidal(NPLMAX)
      real*8 rtidal(NPLMAX), j2tidal(NPLMAX)
      real*8 mu(NPLMAX), eta(NPLMAX), dt, mm, mc, mp
      integer oloc(NPLMAX,NPLMAX)
      integer nbod
      real*8 xj(NPLMAX), yj(NPLMAX), zj(NPLMAX)

c...  Inputs and Outputs
      real*8 sx(NPLMAX), sy(NPLMAX), sz(NPLMAX), s

c...  Internals:
      real*8 gm, r2, r, j2, rs
      real*8 sxnew, synew, sznew
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

          sxnew = sx(itidal) + dt*3d0*gm/(r**4)*j2/(atidal(itidal)*s)*
     &          rs*(yj(j)*sz(itidal)-zj(j)*sy(itidal))
          synew = sy(itidal) + dt*3d0*gm/(r**4)*j2/(atidal(itidal)*s)*
     &          rs*(-xj(j)*sz(itidal)+zj(j)*sx(itidal))
          sznew = sz(itidal) + dt*3d0*gm/(r**4)*j2/(atidal(itidal)*s)*
     &          rs*(xj(j)*sy(itidal)-yj(j)*sx(itidal))

          sx(itidal) = sxnew
          sy(itidal) = synew
          sz(itidal) = sznew

        end if
      end do

      return
      end    ! tidal_quadrupole_spin_evolve
c-----------------------------------------------------------------------
