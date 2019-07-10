c***********************************************************************
c                            CE_CHANGE_ORBIT_TP.F
c***********************************************************************
c Check if a tp should change orbit
c
c Inputs:
c    nbod          ==> Number of massive bodies
c    mass          ==> Masses of bodies
c    oloc          ==> Link Bodies <--> Orbits
c    eta, mu       ==> Masses for centers & sats
c    xb, yb, zb    ==> Bodies position in Barycentric coord
c    xbt, ybt, zbt ==> Tp position in Barycentric coord
c
c Outputs:
c    orbct         ==> Centers of tp's new orbit

      subroutine ce_change_orbit_tp(nbod, mass, eta, mu, oloc, orbct,
     &   xb, yb, zb, xbt, ybt, zbt)

      include '../swift.inc'

c...  Inputs Only:
      integer nbod, oloc(NPLMAX, NPLMAX)
      real*8 mass(nbod), eta(nbod), mu(nbod)
      real*8 xb(nbod), yb(nbod), zb(nbod)
      real*8 xbt, ybt, zbt

c...  Output Only:
      integer orbct(nbod)

c...  Internals:
      integer i, j, jmax
      real*8 xg, yg, zg
      real*8 a(2*nbod-1), r2, maxa

c----
c...  Executable code

      orbct = 0

      do i=1,nbod
         r2 = (xb(i)-xbt)*(xb(i)-xbt)+(yb(i)-ybt)*(yb(i)-ybt)
     &        +(zb(i)-zbt)*(zb(i)-zbt)
         a(i) = mass(i)/r2
      end do

      do j=2,nbod

         xg = 0.d0
         yg = 0.d0
         zg = 0.d0
         do i=1,nbod
            if (oloc(j,i).ne.0) then
               xg = xg + mass(i)*xb(i)/(eta(j)+mu(j))
               yg = yg + mass(i)*yb(i)/(eta(j)+mu(j))
               zg = zg + mass(i)*zb(i)/(eta(j)+mu(j))
            end if
         end do
         r2 = (xg-xbt)*(xg-xbt)+(yg-ybt)*(yg-ybt)+(zg-zbt)*(zg-zbt)
         a(nbod+j-1) = (eta(j)+mu(j))/r2

      end do

      jmax = -1
      maxa = -1.
      do j=1,2*nbod-1
         if (a(j).gt.maxa) then
            jmax = j
            maxa = a(j)
         end if
      end do

      if (jmax.le.nbod) then
         orbct(jmax) = -1
      else
         jmax = jmax-nbod+1
         do i=1,nbod
            if (oloc(jmax,i).ne.0.) orbct(i) = -1
         end do
      end if

      end ! ce_change_orbit_tp
c-----------------------------------------------------------------------
