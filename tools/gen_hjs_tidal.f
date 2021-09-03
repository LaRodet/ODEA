c***********************************************************************
c                         GEN_HJS_TIDAL.F
c***********************************************************************
c This code generates the input file for ODEA from the input of the orbital elements

      program gen_hjs

      include '../swift.inc'

      real*8, parameter :: smassyr = twopi*twopi,
     &                       dr = 1.7453292519943294d-2,      ! pi/180
     &                       au = 1.495978707d8  ,             ! km
     &                       yr = 365.25d0                       ! day

      real*8   mass(nplmax), mtot, rpl(nplmax), rplsq(nplmax)
      real*8   atidal(nplmax), qtidal(nplmax), j2tidal(nplmax)
      real*8   rtidal(nplmax), spinperiod, spin
      real*8   stidalx(nplmax), stidaly(nplmax), stidalz(nplmax)
      real*8   mat(nplmax,nplmax), umat(nplmax,nplmax)
      real*8   eta(nplmax), mu(nplmax), vsat, vcen
      real*8   xb(nplmax), yb(nplmax), zb(nplmax)
      real*8   vxb(nplmax), vyb(nplmax), vzb(nplmax)
      real*8   xj(nplmax), yj(nplmax), zj(nplmax)
      real*8   vxj(nplmax), vyj(nplmax), vzj(nplmax)
      real*8   gm, a, e, inc, capom, omega, capm, q
      real*8   fachill, r2hill(nplmax)

      real*8   etatpt,matpt(nplmax),umatpt(nplmax),masstemp(nplmax)
      real*8   xjtemp(nplmax),yjtemp(nplmax),zjtemp(nplmax)
      real*8   vxjtemp(nplmax),vyjtemp(nplmax),vzjtemp(nplmax)
      real*8   xbtemp(nplmax),ybtemp(nplmax),zbtemp(nplmax)
      real*8   vxbtemp(nplmax),vybtemp(nplmax),vzbtemp(nplmax)
      real*8   umpart(nplmax,nplmax),matp(nplmax,ntpmax)
      real*8   umatp(nplmax,ntpmax)

      real*8   xjt(ntpmax),yjt(ntpmax),zjt(ntpmax)
      real*8   vxjt(ntpmax),vyjt(ntpmax),vzjt(ntpmax)
      real*8   matr(3,3),mati(3,3),idisk,odisk
      real*8   vxh(nplmax),vyh(nplmax),vzh(nplmax)
      real*8   rstat(ntpmax,nstatr)
      real*8   cimax,imax,emax,emin,xx,aminx,amaxx,amin,amax,s
      real*8   xt,yt,zt,vxt,vyt,vzt

      logical lclose,ok,okc,sat,cen

      real*4 seed

      integer   istat(ntpmax,nstat),oloc(nplmax,nplmax),iflgchk
      integer   orbct(nplmax),nc,kk,jj,oloctt(nplmax)
      integer   oloct(nplmax,ntpmax),ntp,ntpt
      integer   ii,k,iseed,nbod,ialpha,i,j
      integer   icflg,irflg,iuflg,dkloc

      character*1    rep
      character*32   name
      character*80   inplfile,intpfile,tidalfile

      ok=.false.
      do while(.not.ok)
         write(*,*) ' units menu (g=1) : '
         write(*,*) '       0 ==> solar masses and au '
         write(*,*) '       1 ==> au and years '
         read(*,*) iuflg
         ok = ((iuflg.eq.0).or.(iuflg.eq.1))
      end do

      ok=.false.
      do while(.not.ok)
         write(*,*) ' coordinate menu: '
         write(*,*) '       0 ==> ecliptic '
         write(*,*) '       1 ==> invariable plane '
         read(*,*) icflg
         ok = ((icflg.eq.0).or.(icflg.eq.1))
      end do

      ok=.false.
      do while(.not.ok)
         write(*,*) ' planetary radius menu: '
         write(*,*) '  0 ==> do not include a radius '
         write(*,*) '  1 ==> include the physical radius of planet '
         write(*,*) '  2 ==> include a multiple of the hill radius '
         read(*,*) irflg
         ok = ((irflg.eq.0).or.(irflg.eq.1).or.(irflg.ne.2))
      end do

      if(irflg.eq.2) then
         write(*,*) ' what is that multiple ? '
         read(*,*) fachill
      end if

      write(*,*)' give number of massive bodies in the system '
      read(*,*) nbod
      do i=1,nbod
          write(*,*)" give mass of body #",i," in solar masses and ",
     & "tidal parameters alpha, Q', radius (km), spin period (d), ",
     & "spin orientation sx sy sz and J2"
          read(*,*) mass(i), atidal(i), qtidal(i), rtidal(i),
     &        spinperiod, stidalx(i), stidaly(i), stidalz(i), j2tidal(i)
          rtidal(i) = rtidal(i)/au
          if (spinperiod.gt.0d0) then
            spin = (stidalx(i)**2 + stidaly(i)**2 + stidalz(i)**2)
            if (spin.ne.1d0) then
              write(*,*)"Fatal error: Wrong spin normalization"
              stop
            end if
          spin = 2d0*PI/(spinperiod/yr)
          stidalx(i) = spin*stidalx(i)
          stidaly(i) = spin*stidaly(i)
          stidalz(i) = spin*stidalz(i)
        else
          stidalx(i) = 0d0
          stidaly(i) = 0d0
          stidalz(i) = 0d0
        end if
      end do
      if (irflg.eq.1) then
         write(*,*)' give radius of body #',i,' in km'
         do i=1,nbod
            read(*,*) rpl(i)
         end do
      end if

      if(iuflg.eq.1) mass(1:nbod) = mass(1:nbod)*smassyr

      xj(1) = 0.0
      yj(1) = 0.0
      zj(1) = 0.0
      vxj(1) = 0.0
      vyj(1) = 0.0
      vzj(1) = 0.0

      do i=2,nbod
         write(*,*)' Information relative to orbit #',i-1,':'
         write(*,*)' situation of bodies relative to that orbit?'
         write(*,*)' (0 = foreign,   1 = satellite,  -1 = cente)r'
         read(*,*)oloc(i,1:nbod)

         write(*,*) ' orbital parameters?'
         read(*,*)a, e, inc, capom, omega, capm

         capm = dr*capm
         omega = dr*omega
         capom = dr*capom
         inc = dr*inc

         gm = 0.
         do j = 1,nbod
            if (oloc(i,j).ne.0) gm = gm+mass(j)
         end do
         if (e.lt.1) then
            ialpha = -1
         else if (e.eq.1) then
            ialpha = 0
         else
            ialpha = 1
         end if
         call orbel_el2xv(gm, ialpha, a, e, inc, capom, omega, capm,
     &      xj(i), yj(i), zj(i), vxj(i), vyj(i), vzj(i))
      end do

c...  Calc the radii
      if(irflg.eq.0) then
         lclose = .false.
         rplsq(1:nbod) = 0.d0
      else if(irflg.eq.1) then
         lclose = .true.
         do i=1,nbod
            rplsq(i) = (rpl(i)/au)**2
         end do
      endif

c... Compute eta's and mu's: center and satellite masses for orbits
      eta = 0
      mu = 0
      do j = 2,nbod
         do i = 1,nbod
            if (oloc(j,i).eq.1) mu(j) = mu(j)+mass(i)
            if (oloc(j,i).eq.-1) eta(j) = eta(j)+mass(i)
         end do
      end do

c... Build transform matrix barycentric --> jacobi
      mat = 0.0d0
      mtot = sum(mass(1:nbod))
      do i = 1,nbod
         mat(1,i) = mass(i)/mtot
      end do
      do j = 2,nbod
         do i = 1,nbod
            if (oloc(j,i).eq.1) mat(j,i) = mass(i)/mu(j)
            if (oloc(j,i).eq.-1) mat(j,i) = -mass(i)/eta(j)
         end do
      end do

c... Build inverse transform matrix jacobi --> barycentric
      umat = 0.0d0
      do i = 1,nbod
         umat(i,1) = 1.0d0
      end do
      do j = 2,nbod
         vsat = eta(j)/(mu(j)+eta(j))
         vcen = -mu(j)/(mu(j)+eta(j))
         do i = 1,nbod
            if (oloc(j,i).eq.1) umat(i,j) = vsat
            if (oloc(j,i).eq.-1) umat(i,j) = vcen
         end do
      end do

      write(*,*) ' '
      write(*,*) ' Name of planet data file?'
      read(*,'(a)') inplfile
      iflgchk = 0
      call io_dump_pl_hjs(inplfile, nbod, oloc, mass(1:nbod), umat,
     &   xj(1:nbod), yj(1:nbod), zj(1:nbod), vxj(1:nbod), vyj(1:nbod),
     &   vzj(1:nbod), lclose, iflgchk, rplsq(1:nbod))
      write(*,*) ' Name of tidal data file?'
      read(*,'(a)') tidalfile
      call io_dump_tidal(tidalfile, nbod, mass(1:nbod), atidal(1:nbod),
     &   qtidal(1:nbod), rtidal(1:nbod), stidalx(1:nbod),
     &   stidaly(1:nbod), stidalz(1:nbod), j2tidal(1:nbod))

      okc = .true.
      write(*,*) ' Input iseed (large and odd):'
      read(*,*) iseed
      seed=float(iseed)
      call random_seed(iseed)
      ntpt = 0
      ii = 0
      do while(okc)
         write(*,*) ' Number of test particles?'
         read(*,*) ntp
         okc = (ntp.ne.0)
         if (okc) then
            ntpt = ntpt+ntp
            write(*,*) 'Centers for tp?'
            read(*,*) (orbct(i),i=1,nbod)
            write(*,*) 'Centers: ',(orbct(j),j=1,nbod)
            write(*,*) 'Disk in ecliptic (0), centers (1)',
     &         ', invariable (2), other (3), or orbit (<0) ?'
            read(*,*) dkloc

c... Check the validity of the hierarchy
            call hierarchktp(nbod,oloc,orbct,ok)
            write(*,*)' '
            if (.not.ok) then
               write(*,*)'Stopping: bad hiearchical configuration'
               stop
            end if
            write(*,*)'Hierarchical structure ok!'
            write(*,*)' '

c... Computation of the oloct array (hierarchy of the disk)
            do j=1,nbod
               oloctt(j) = 0
            end do
            do j = 2,nbod
               sat = .true.
               cen = .true.
               do i = 1,nbod
                  if (orbct(i).eq.-1) then
                     sat = sat.and.(oloc(j,i).eq.1)
                     cen = cen.and.(oloc(j,i).eq.-1)
                  end if
               end do
               if (sat) oloctt(j) = 1
               if (cen) oloctt(j) = -1
            end do
            print*,'oloct',(oloctt(j),j=1,nbod)

            etatpt = 0.0d0
            matpt(:) = 0.0d0
            umatpt(:) = 0.0d0
            do j = 1,nbod
               if (orbct(j).eq.-1) etatpt = etatpt+mass(j)
            end do
            do j = 1,nbod
               if (orbct(j).eq.-1) matpt(j) = -mass(j)/etatpt
            end do
            do k = 1,nbod
               if (orbct(k).eq.-1) then
                  do j=1,nbod
                     umatpt(j) = umatpt(j) - matpt(k)*umat(k,j)
                  end do
               end if
            end do

            write(*,*) ' Eccentricities of particles (min, max)?'
            read(*,*) emin, emax
            write(*,*) ' Max inclination of particles (in deg)?'
            read(*,*) imax
            write(*,*) ' Semi-major axes (min, max)?'
            read(*,*) amin, amax
            write(*,*) ' Power index?'
            read(*,*) xx
            aminx = amin**(xx+1.)
            amaxx = amax**(xx+1.)
            e = emin
            cimax = cos(imax/degrad)

            nc = 0
            do j = 1,nbod
               if (orbct(j).eq.-1) then
                  nc = nc+1
               end if
            end do

c... Computation of the tranformation matrix to the tp's plane
c...    Tps in ecliptic plane
            if ((dkloc.eq.0).or.((dkloc.eq.1).and.(nc.eq.1))) then

               call invar(1, [1.d0], [1.d0], [0.d0], [0.d0], [0.d0],
     &         [1.d0], [0.d0], matr, 0)

c...    Tps in the centers invariable plane (computed only if more than one center)
            elseif ((nc.gt.1).and.(dkloc.eq.1)) then

c...       extraction of a submatrix
               umpart = 0.0d0
               do j = 1,nc
                  umpart(j,1) = 1.0d0
               end do
               jj = 1
               do j = 2,nbod
                  ok = .true.
                  do k = 1,nbod
                     if ((matpt(k).eq.0.0d0).and.(oloc(j,k).ne.0))
     &                  ok = .false.
                  end do
                  if (ok) then
                     jj = jj+1
                     kk = 0
                     xjtemp(jj) = xj(j)
                     yjtemp(jj) = yj(j)
                     zjtemp(jj) = zj(j)
                     vxjtemp(jj) = vxj(j)
                     vyjtemp(jj) = vyj(j)
                     vzjtemp(jj) = vzj(j)
                     do k = 1,nbod
                        if (matpt(k).ne.0.0d0) then
                           kk = kk+1
                           masstemp(kk) = mass(k)
                           umpart(kk,jj) = umat(k,j)
                        end if
                     end do
                  end if
               end do

               call coord_g2b(nc, umpart, masstemp(1:nc), xjtemp(1:nc),
     &            yjtemp(1:nc), zjtemp(1:nc), vxjtemp(1:nc),
     &            vyjtemp(1:nc), vzjtemp(1:nc), xbtemp(1:nc),
     &            ybtemp(1:nc), zbtemp(1:nc), vxbtemp(1:nc),
     &            vybtemp(1:nc), vzbtemp(1:nc))

               call invar(nc, masstemp(1:nc), xbtemp(1:nc),
     &            ybtemp(1:nc), zbtemp(1:nc), vxbtemp(1:nc),
     &            vybtemp(1:nc), vzbtemp(1:nc), matr, 1)

c...    Tps in the total invariable plan
            elseif ((dkloc.eq.2)) then

               call coord_g2b(nbod, umat, mass(1:nbod), xj(1:nbod),
     &         yj(1:nbod), zj(1:nbod), vxj(1:nbod), vyj(1:nbod),
     &         vzj(1:nbod), xb(1:nbod), yb(1:nbod), zb(1:nbod),
     &         vxb(1:nbod), vyb(1:nbod), vzb(1:nbod))

               call invar(nbod,mass,xb,yb,zb,vxb,vyb,vzb,matr,0)

c...    Tps in the plane of one given orbit
            elseif (dkloc.lt.0) then

               call invar(1, [1.0d0], xj(-dkloc), yj(-dkloc),
     &          zj(-dkloc), vxj(-dkloc), vyj(-dkloc), vzj(-dkloc), matr,
     &          0)

c...    Tps in a plane defined by two angles
            elseif (dkloc.eq.3) then

               write(*,*) ' Disk inclination (in deg)?'
               read(*,*) idisk
               write(*,*) ' Disk position angle (in deg)?'
               read(*,*) odisk
               idisk = idisk*dr
               odisk = odisk*dr

               matr(1,1) = cos(odisk)
               matr(1,2) = sin(odisk)
               matr(1,3) = 0.0d0

               matr(2,1) = -cos(idisk)*sin(odisk)
               matr(2,2) = cos(idisk)*cos(odisk)
               matr(2,3) = sin(idisk)

               matr(3,1) = sin(idisk)*sin(odisk)
               matr(3,2) = -sin(idisk)*cos(odisk)
               matr(3,3) = cos(idisk)

               open(17, file='matpass.dat', status='unknown')
               write(17,*) matr(1,1), matr(1,2), matr(1,3)
               write(17,*) matr(2,1), matr(2,2), matr(2,3)
               write(17,*) matr(3,1), matr(3,2), matr(3,3)
               close(17)

            end if

c... Generation of the Tps

            write(*,*) ' Generation of the test particles...'
            do i=1,ntp
               ii = ii+1
               ialpha = -1

               call random_number(s)
               inc = acos( 1.0d0 - (1.0d0-cimax)*s)
               call random_number(s)
               capm = twopi*s
               call random_number(s)
               omega = twopi*s
               call random_number(s)
               capom = twopi*s
               call random_number(s)
               a = (aminx + (amaxx - aminx)*s)**(1./(xx+1.))
               call random_number(s)
               e = emin + (emax - emin)*s

               call orbel_el2xv(etatpt, ialpha, a, e, inc, capom, omega,
     &            capm, xjt(ii), yjt(ii), zjt(ii), vxjt(ii), vyjt(ii),
     &            vzjt(ii))

               xt = matr(1,1)*xjt(ii)+matr(2,1)*yjt(ii)
     &            +matr(3,1)*zjt(ii)
               yt = matr(1,2)*xjt(ii)+matr(2,2)*yjt(ii)
     &            +matr(3,2)*zjt(ii)
               zt = matr(1,3)*xjt(ii)+matr(2,3)*yjt(ii)
     &            +matr(3,3)*zjt(ii)
               xjt(ii) = xt
               yjt(ii) = yt
               zjt(ii) = zt
               vxt = matr(1,1)*vxjt(ii)+matr(2,1)*vyjt(ii)
     &            +matr(3,1)*vzjt(ii)
               vyt = matr(1,2)*vxjt(ii)+matr(2,2)*vyjt(ii)
     &            +matr(3,2)*vzjt(ii)
               vzt = matr(1,3)*vxjt(ii)+matr(2,3)*vyjt(ii)
     &            +matr(3,3)*vzjt(ii)
               vxjt(ii) = vxt
               vyjt(ii) = vyt
               vzjt(ii) = vzt

               oloct(1:nbod,ii) = oloctt(1:nbod)
               matp(1:nbod,ii) = matpt(1:nbod)
               umatp(1:nbod,ii) = umatpt(1:nbod)

               do j=1,nstat
                  istat(ii,j) = 0
               end do
               do j=1,nstatr
                  rstat(ii,j) = 0.0d0
               end do
            end do
         end if
      end do

      write(*,*) 'Name of test particle data file?'
      read(*,'(a)') intpfile
      print*, intpfile

      call io_dump_tp_hjs(intpfile, nbod, ntpt, matp,
     &   xjt(1:ntpt), yjt(1:ntpt), zjt(1:ntpt), vxjt(1:ntpt),
     &   vyjt(1:ntpt), vzjt(1:ntpt), istat, rstat)

      end program gen_hjs

c***********************************************************************
c                         Invar
c***********************************************************************
c Computation of the transformation matrix to the invariable plane
      subroutine invar(nbod, mass, x, y, z, vx, vy, vz, a, flag)

      include '../swift.inc'

      integer i, nbod

      real*8 mass(nbod), x(nbod), y(nbod), z(nbod)
      real*8 vx(nbod), vy(nbod), vz(nbod)
      real*8 xt, yt, zt, vxt, vyt, vzt
      real*8 c1, c2, c3, c, bot, a(3,3)

      integer flag ! 1 to compute x,y,z,vx,vy,vz in the new plane


      c1 = sum(mass(1:nbod)*
     &   (y(1:nbod)*vz(1:nbod)-z(1:nbod)*vy(1:nbod)))
      c2 = sum(mass(1:nbod)*
     &   (z(1:nbod)*vx(1:nbod)-x(1:nbod)*vz(1:nbod)))
      c3 = sum(mass(1:nbod)*
     &   (x(1:nbod)*vy(1:nbod)-y(1:nbod)*vx(1:nbod)))

      c = sqrt( c1*c1 + c2*c2 + c3*c3 )
      c1 = c1/c
      c2 = c2/c
      c3 = c3/c

      print*, 'Normal vector to invariable plane', c1, c2, c3

      bot = 1.0d0/(1.0d0 + c3)

      a(1,1) = 1.0d0 - c1*c1*bot
      a(1,2) = -1.0d0*c1*c2*bot
      a(1,3) = -1.0d0*c1

      a(2,1) = -1.0d0*c1*c2*bot
      a(2,2) = 1.0d0 - c2*c2*bot
      a(2,3) = -1.0d0*c2

      a(3,1) = c1
      a(3,2) = c2
      a(3,3) = c3

      if (flag.eq.1) then

         do i=1,nbod
            xt = a(1,1)*x(i) + a(1,2)*y(i) + a(1,3)*z(i)
            yt = a(2,1)*x(i) + a(2,2)*y(i) + a(2,3)*z(i)
            zt = a(3,1)*x(i) + a(3,2)*y(i) + a(3,3)*z(i)
            x(i) = xt
            y(i) = yt
            z(i) = zt
            vxt = a(1,1)*vx(i) + a(1,2)*vy(i) + a(1,3)*vz(i)
            vyt = a(2,1)*vx(i) + a(2,2)*vy(i) + a(2,3)*vz(i)
            vzt = a(3,1)*vx(i) + a(3,2)*vy(i) + a(3,3)*vz(i)
            vx(i) = vxt
            vy(i) = vyt
            vz(i) = vzt
         end do

      end if

      open(17, file='matpass.dat', status='unknown')
      write(17,*) a(1,1), a(1,2), a(1,3)
      write(17,*) a(2,1), a(2,2), a(2,3)
      write(17,*) a(3,1), a(3,2), a(3,3)
      close(17)

      end ! invar

c***********************************************************************
c        Hierarchktp
c***********************************************************************
c Subroutine dedicated to verify the hierarchical structure of the tp's

      subroutine hierarchktp(nbod, oloc, orbct, ok)

      include '../swift.inc'

      integer i,k,l,nnn,nbod,oloc(nplmax,nplmax),orbct(nplmax)
      logical test,ok

      write(*,*)' '
      ok = .true.

c... First check that the orbct array is made only of 0's and -1's
      do i=1,nbod
         ok = ok.and.((orbct(i).eq.0).or.(orbct(i).eq.-1))
      end do
      if (.not.ok) then
         write(*,*)'the orbit given is not a valid tp orbit'
         return
      end if

c... Adding a tp means adding an orbit. we check the validity of the
c... the new orbit in the hierarchical structure of the system, i.e.
c... with all orbits i within massive bodies
      do i = 2,nbod
         test = .true.
c... We first test whether tp orbit and orbit i are foreign
         do k = 1,nbod
            test = test.and.(oloc(i,k)*orbct(k).eq.0)
         end do
         if (.not.test) then
c... From here on they are not foreign. We check whether tp orbit
c... is inner to orbit i, i.e. all centers of tp orbit
c... must fall in the centers or the satellites of orbit i.
            l = 1
            do while(orbct(l).eq.0)
               l = l+1
            end do
            nnn = oloc(i,l)
            if (nnn.ne.0) then
               test = .true.
               do k = l+1,nbod
                  if (orbct(k).eq.-1) test =
     &               test.and.(oloc(i,k).eq.nnn)
               end do
            end if
            if ((nnn.eq.0).or.(.not.test)) then
c... Now tp orbit is not inner to orbit i. We check whether orbit i
c... is inner to tp orbit
               l = 1
               do while(oloc(i,l).eq.0)
                  l = l+1
               end do
               nnn = orbct(l)
               if (nnn.ne.0) then
                  test = .true.
                  do k = l+1,nbod
                     if (oloc(i,k).ne.0) test =
     &                  test.and.(orbct(k).eq.nnn)
                  end do
                  if (.not.test) then
c... If all fails then we are sure that orbits i and j do not
c... fulfill the hierarchy rules
                     write(*,*)'Structure problem relative to orbit',i-1
                     ok = .false.
                  end if
               end if
            end if
         end if
      end do

      end ! hierarchktp
