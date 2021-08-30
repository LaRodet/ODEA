c***********************************************************************
c                           ORBEL_XV2EL.F
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

      subroutine orbel_xv2el(x, y, z, vx, vy, vz, gmsum, ialpha,
     &   a, e, inc, capom, omega, capm)

      include '../swift.inc'

c...  Inputs Only:
      real*8 x, y, z, vx, vy, vz, gmsum

c...  Outputs
      integer ialpha
      real*8 a, e, inc, capom, omega, capm

c...  Internals:
      real*8 hx, hy, hz, h2, h, r, v2, v, vdotr, energy
      real*8 fac, face, cape, capf, tmpf
      real*8 cw, sw, w, u

c----
c...  Executable code

c... Compute the angular momentum H, and thereby the inclination INC.

      hx = y*vz - z*vy
      hy = z*vx - x*vz
      hz = x*vy - y*vx
      h2 = hx*hx + hy*hy +hz*hz
      h  = sqrt(h2)
      inc = acos(hz/h)

c...  Compute longitude of ascending node CAPOM
c     and the argument of latitude u.
      fac = (hx**2 + hy**2)/(h**2)

      if(fac.lt. TINY ) then
         capom = 0.d0
         u = atan2(y,x)
         if(abs(inc - PI).lt. 10.d0*TINY) u = -u
      else
         capom = atan2(hx,-hy)
         u = atan2 ( z/sin(inc) , x*cos(capom) + y*sin(capom))
      endif

      if(capom .lt. 0.d0) capom = capom + 2.d0*PI
      if(u .lt. 0.d0) u = u + 2.d0*PI

c...  Compute the radius R and velocity squared V2, the dot
c...  product RDOTV, and the energy per unit mass ENERGY.

      r = sqrt(x*x + y*y + z*z)
      v2 = vx*vx + vy*vy + vz*vz
      v = sqrt(v2)
      vdotr = x*vx + y*vy + z*vz
      energy = 0.5d0*v2 - gmsum/r

c...  Determine type of conic section and label it via IALPHA
      if(abs(energy*r/gmsum) .lt. sqrt(TINY)) then
         ialpha = 0
      else
         if(energy .lt. 0.d0) ialpha = -1
         if(energy .gt. 0.d0) ialpha = +1
      endif

c... Depending on the conic type, determine the remaining elements

c... ELLIPSE :
      if(ialpha .eq. -1) then
         a = -0.5d0*gmsum/energy
         fac = 1.d0 - h2/(gmsum*a)

         if (fac .gt. TINY) then
            e = sqrt ( fac )
            face =(a-r)/(a*e)

            if ( face .gt. 1.d0) then
               cape = 0.d0
            else
               if ( face .gt. -1.d0) then
                  cape = acos( face )
               else
                  cape = PI
               endif
            endif

            if ( vdotr .lt. 0.d0 ) cape = 2.d0*PI - cape
            cw = (cos( cape) -e)/(1.d0 - e*cos(cape))
            sw = sqrt(1.d0 - e*e)*sin(cape)/(1.d0 - e*cos(cape))
            w = atan2(sw,cw)
            if(w .lt. 0.d0) w = w + 2.d0*PI
         else
            e = 0.d0
            w = u
            cape = u
         endif

         capm = cape - e*sin (cape)
         omega = u - w
         if(omega .lt. 0.d0) omega = omega + 2.d0*PI
         omega = omega - int(omega/(2.d0*PI))*2.d0*PI

      endif

c... HYPERBOLA
      if(ialpha .eq. +1) then

         a = +0.5d0*gmsum/energy
         fac = h2/(gmsum*a)

         if (fac .gt. TINY) then
            e = sqrt ( 1.d0 + fac )
            tmpf = (a+r)/(a*e)
            if(tmpf.lt.1.0d0) then
               tmpf = 1.0d0
            endif
            capf = log(tmpf + sqrt(tmpf*tmpf -1.d0))
            if ( vdotr .lt. 0.d0 ) capf = - capf
            cw = (e - cosh(capf))/(e*cosh(capf) - 1.d0 )
            sw = sqrt(e*e - 1.d0)*sinh(capf)/(e*cosh(capf) - 1.d0 )
            w = atan2(sw,cw)
            if(w .lt. 0.d0) w = w + 2.d0*PI
         else
c... We only get here if a hyperbola is essentially a parabola
c...  so we calculate e and w accordingly to avoid singularities
            e = 1.d0
            tmpf = 0.5d0*h2/gmsum
            w = acos(2.d0*tmpf/r -1.d0)
            if ( vdotr .lt. 0.d0) w = 2.d0*PI - w
            tmpf = (a+r)/(a*e)
            capf = log(tmpf + sqrt(tmpf*tmpf -1.d0))
         endif

         capm = e * sinh(capf) - capf
         omega = u - w
         if(omega .lt. 0.d0) omega = omega + 2.d0*PI
         omega = omega - int(omega/(2.d0*PI))*2.d0*PI
      endif

c... PARABOLA (in this case we use "a" to mean pericentric distance)
      if(ialpha .eq. 0) then
         a =  0.5d0*h2/gmsum
         e = 1.d0
         w = acos(2.d0*a/r -1.d0)
         if ( vdotr .lt. 0.d0) w = 2.d0*PI - w
         tmpf = tan(0.5d0 * w)
         capm = tmpf* (1.d0 + tmpf*tmpf/3.d0)
         omega = u - w
         if(omega .lt. 0.d0) omega = omega + 2.d0*PI
         omega = omega - int(omega/(2.d0*PI))*2.d0*PI
      endif

      return
      end    ! orbel_xv2el
c-----------------------------------------------------------------------
