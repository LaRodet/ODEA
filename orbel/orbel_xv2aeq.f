c***********************************************************************
c                            ORBEL_XV2AEQ.F
c***********************************************************************
c Purpose:  Given the cartesian position and velocity of an orbit,
c       compute the osculating orbital elements a, e, and q only.
c
c Inputs:
c    x, y, z    ==>  position of object
c    vx, vy, vz ==>  velocity of object
c    gmsum      ==> G * total mass
c
c Outputs:
c    ialpha   ==> conic section type
c    a        ==> semi-major axis or pericentric distance if a parabola
c    e        ==> eccentricity
c    q        ==> perihelion distance
c
c Algorithm: See e.g. p.70 of Fitzpatrick's "Priciples of Cel. Mech."
c Remarks: Based on M. Duncan's orbel_xv2el.f
c      This routine is generally applied to study (hyperbolic) close
c      encounters of test particles with planets.

      subroutine orbel_xv2aeq(x, y, z, vx, vy, vz, gmsum,
     &   ialpha, a, e, q)

      include '../swift.inc'

c...  Inputs Only:
      real*8 x, y, z, vx, vy, vz, gmsum

c...  Outputs
      integer ialpha
      real*8 a, e, q

c...  Internals:
      real*8 hx, hy, hz, h2, r, v2, energy, fac
c----
c...  Executable code

c...  Compute the angular momentum H, and thereby the inclination INC.
      hx = y*vz - z*vy
      hy = z*vx - x*vz
      hz = x*vy - y*vx
      h2 = hx*hx + hy*hy + hz*hz

c...  Compute the radius R and velocity squared V2, and the dot
c...  product RDOTV, the energy per unit mass ENERGY .

      r = sqrt(x*x + y*y + z*z)
      v2 = vx*vx + vy*vy + vz*vz
      energy = 0.5d0*v2 - gmsum/r

c...  Determine type of conic section and label it via IALPHA
      if(abs(energy*r/gmsum) .lt. sqrt(TINY)) then
         ialpha = 0
      else
         if(energy .lt. 0.d0) ialpha = -1
         if(energy .gt. 0.d0) ialpha = +1
      endif

c...  Depending on the conic type, determine the remaining elements

c...  Ellipse:
      if(ialpha .eq. -1) then
         a = -0.5d0*gmsum/energy
         fac = 1.d0 - h2/(gmsum*a)

         if (fac .gt. TINY) then
            e = sqrt ( fac )
         else
            e = 0.d0
         endif

         q = a*(1.d0 - e)
      endif

c...  Hyperbola
      if(ialpha .eq. +1) then
         a = +0.5d0*gmsum/energy
         fac = h2/(gmsum*a)
         if (fac .gt. TINY) then
            e = sqrt ( 1.d0 + fac )
            q = -a*(1.d0 - e)
c...     Have to insert minus sign in expression for q because this code
c...     takes a > 0, even for a hyperbola
         else
c...     We only get here if a hyperbola is essentially a parabola
c...     so we calculate e accordingly to avoid singularities
            e = 1.d0
            q = 0.5*h2/gmsum
         endif
      endif

c...  Parabola: ( NOTE - in this case "a", which is formally infinite,
c         is arbitrarily set equal to the pericentric distance q).
      if(ialpha .eq. 0) then
         a =  0.5d0*h2/gmsum
         e = 1.d0
         q = a
      endif

      return
      end    ! orbel_xv2aeq
c-----------------------------------------------------------------------
