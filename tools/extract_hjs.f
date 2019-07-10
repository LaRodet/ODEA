c***********************************************************************
c                             EXTRACT_HJS.F
c***********************************************************************
c Compute the coordinates of bodies and tp from the bin file
c From a bin file with the orbital elements, create 4 files:
c    xvbodies.dat  ==> Barycentric positions massive bodies
c    xvtp.dat      ==> Barycentric positions tp
c    elbodies.dat  ==> Orbital elements massive bodies
c    eltp.dat      ==> Orbital elements tp

      program extract_hjs

      include '../swift.inc'

      real*8,parameter :: b = 0d0 ! blank

      real*8 t, t0, tstop, dt, dtout, dtdump, dtpout, ttp
      real*8 verbose, tverbose
      real*8 t0tp, tstoptp
      real*8 rmin, rmax, rmaxu, qmin, rplsq
      real*8 a, e, inc, capom, omega, capm, gm
      real*8 mat(NPLMAX,NPLMAX), umat(NPLMAX,NPLMAX)
      real*8 mass(NPLMAX), eta(NPLMAX), mu(NPLMAX)
      real*8 matp(NPLMAX,NTPMAX), umatp(NPLMAX,NTPMAX)
      real*8 etatp(NTPMAX)

      real*8 xj(NPLMAX), yj(NPLMAX), zj(NPLMAX)
      real*8 vxj(NPLMAX), vyj(NPLMAX), vzj(NPLMAX)
      real*8 xb(NPLMAX), yb(NPLMAX), zb(NPLMAX)
      real*8 vxb(NPLMAX), vyb(NPLMAX), vzb(NPLMAX)

      real*8 xjt(NTPMAX), yjt(NTPMAX), zjt(NTPMAX)
      real*8 vxjt(NTPMAX), vyjt(NTPMAX), vzjt(NTPMAX)
      real*8 xt, yt, zt, vxt, vyt, vzt
      real*8 xbt, ybt, zbt
      real*8 vxbt, vybt, vzbt

      real*8 matr(3,3)

      integer iflgchk, iu, iread, ireadtp, i, j, k
      integer nbod, nleft, id, ntp
      integer io_read_hdr_odea, io_read_line, ierr
      integer io_read_line_r
      integer oloc(NPLMAX,NPLMAX), ialpha
      integer oloct(NPLMAX,NTPMAX)
      integer istat(NTPMAX,NSTAT), rstat(NTPMAX,NSTATR)

      logical lclose

      character*80 inparfile, diro, dirs, gname, outfile, fopenstat
      character*80 inplfile, intpfile, matfile
      character*80 xvtp, eltp, xvbodies, elbodies

 1    format('t =',F15.3,' nbod =',I2)

      write(*,*) 'enter name of parameter data file : '
      read(*,'(a)') inparfile
      call io_init_param(inparfile, t0, tstop, dt, dtout, dtdump,
     &   iflgchk, rmin, rmax, rmaxu, qmin, lclose, diro, dirs,
     &   gname, outfile, fopenstat)

      verbose = (tstop-t0)/10.

      write(*,*) 'enter name of planet data file : '
      read(*,'(a)') inplfile
      call io_init_pl_hjs(inplfile, lclose, nbod, oloc,
     &   mass, eta, mu, mat, umat, xj, yj, zj, vxj, vyj, vzj, rplsq)

      write(*,*) 'Enter name of test particle data file : '
      read(*,'(a)') intpfile
      call io_init_tp_hjs(intpfile, nbod, ntp, oloc, oloct, mass, umat,
     &   etatp, matp, umatp, xjt, yjt, zjt, vxjt, vyjt, vzjt, istat,
     &   rstat)

      iu = 20

      open(unit=iu, file=trim(diro)//'/'//outfile,
     &   status='old',form='unformatted')

      matfile = trim(diro)//'/mat.dat'
      open(0, file=matfile, status='unknown')
      write(0,*) eta(1:nbod)
      write(0,*) mu(1:nbod)
      do j=1,nbod
         write(0,*) mat(j,1:nbod)
      end do
      do j=1,nbod
         write(0,*) umat(j,1:nbod)
      end do
      close(0)

      elbodies = trim(diro)//'/elbodies.dat'
      open(1, file=elbodies, status='unknown')
      write(1,*) nbod

      xvbodies = trim(diro)//'/xvbodies.dat'
      open(2, file=xvbodies, status='unknown')
      write(2,*) nbod

      if (ntp.gt.0) then

         write(*,*) 'Enter time interval for test particles: '
         read(*,*) t0tp, tstoptp, dtpout
         print*, t0tp, tstoptp, dtpout
         ttp = t0tp

         eltp = trim(diro)//'/eltp.dat'
         open(3,file=eltp,status='unknown')
         write(3,*) nbod, ntp

         xvtp = trim(diro)//'/xvtp.dat'
         open(4,file=xvtp,status='unknown')
         write(4,*) nbod, ntp

         open(unit=17, file=trim(diro)//'/'//'matpass.dat',
     &          status='old')
         if(ierr.ne.0) then
            write(*,*) 'Fatal error:  Could not open matpass.dat'
            stop
         endif

         read(17,*) matr(1,1), matr(1,2), matr(1,3)
         read(17,*) matr(2,1), matr(2,2), matr(2,3)
         read(17,*) matr(3,1), matr(3,2), matr(3,3)
         close(17)

      end if

      iread = 0
      ireadtp = 0
      tverbose = 0.
      t = 0.

      do while (t.lt.tstop)

         ierr = io_read_hdr_odea(iu, t, nbod, nleft)

         if (ierr.ne.0) then
            print*,"Problem opening binary file"
            exit
         end if

         iread = iread + 1

         if (t.ge.tverbose) then
            write(6,1)t,nbod
            tverbose = tverbose + verbose
         end if

         write(1,*)t,0.,0.,0.,0.,0.,0.
         write(2,*)t,0.,0.,0.,0.,0.,0.

         do k=2,nbod

            if (btest(iflgchk,0))  then ! bit 0 i set = integer*2
               ierr = io_read_line(iu,id,a,e,inc,capom,omega,capm)
            else
               ierr = io_read_line_r(iu,id,a,e,inc,capom,omega,capm)
            endif

            write(1,*)id,a,e,inc,capom,omega,capm

            gm = eta(k)+mu(k)

            if (e.lt.1) then
               ialpha = -1
            else if (e.eq.1) then
               ialpha = 0
            else
               ialpha = 1
            end if

            call orbel_el2xv(gm, ialpha, a, e, inc, capom, omega, capm,
     &         xj(k), yj(k), zj(k), vxj(k), vyj(k), vzj(k))

         end do

         call coord_g2b(nbod, umat, mass, xj, yj, zj, vxj, vyj, vzj,
     &      xb, yb, zb, vxb, vyb, vzb)

         do i=1,nbod
            write(2,*)i, xb(i), yb(i), zb(i), vxb(i), vyb(i), vzb(i)
         end do

         if ((ntp.gt.0).and.(((t.ge.ttp).and.(t.le.tstoptp))
     &      .or.(t.eq.tstoptp))) then

            ttp = ttp + dtpout

            ireadtp = ireadtp + 1

            write(3,*) t, nleft, 0., 0., 0., 0., 0.
            write(4,*) t, nleft, 0., 0., 0., 0., 0.

            do i=1,nleft
               if(btest(iflgchk,0))  then ! bit 0 is set
                  ierr = io_read_line(iu, id, a, e, inc, capom, omega,
     &               capm)
               else
                  ierr = io_read_line_r(iu, id, a, e, inc, capom, omega,
     &               capm)
               endif

               write(3,*) id, a, e, inc, capom, omega, capm

               gm = etatp(id)

               if (e.lt.1) then
                  ialpha = -1
               else if (e.eq.1) then
                  ialpha = 0
               else
                  ialpha = 1
               end if

               call orbel_el2xv(gm, ialpha, a, e, inc, capom, omega,
     &            capm, xjt(i), yjt(i), zjt(i), vxjt(i), vyjt(i),
     &            vzjt(i))

               xt = matr(1,1)*xjt(i)+matr(2,1)*yjt(i)+matr(3,1)*zjt(i)
               yt = matr(1,2)*xjt(i)+matr(2,2)*yjt(i)+matr(3,2)*zjt(i)
               zt = matr(1,3)*xjt(i)+matr(2,3)*yjt(i)+matr(3,3)*zjt(i)
               vxt = matr(1,1)*vxjt(i)+matr(2,1)*vyjt(i)+matr(3,1)
     &                *vzjt(i)
               vyt = matr(1,2)*vxjt(i)+matr(2,2)*vyjt(i)+matr(3,2)
     &                *vzjt(i)
               vzt = matr(1,3)*vxjt(i)+matr(2,3)*vyjt(i)+matr(3,3)
     &                *vzjt(i)

               call coord_g2b_tp(nbod, umatp(1:nbod,i), xj(1:nbod),
     &            yj(1:nbod), zj(1:nbod), vxj(1:nbod), vyj(1:nbod),
     &            vzj(1:nbod), xt, yt, zt, vxt, vyt, vzt,
     &            xbt, ybt, zbt, vxbt, vybt, vzbt)

               write(4,*) id, xbt, ybt, zbt, vxbt, vybt, vzbt

            end do

         else

            do i=1,nleft
               if(btest(iflgchk,0))  then ! bit 0 is set
                  ierr = io_read_line(iu, id, a, e, inc, capom, omega,
     &               capm)
               else
                  ierr = io_read_line_r(iu, id, a, e, inc, capom, omega,
     &               capm)
               endif
            end do

         end if
      end do

      print*, iread, 'frames'
      print*, ireadtp, 'frames tp'

      close(1)
      close(2)

      if (ntp.gt.0) then
         close(3)
         close(4)
      end if

      end ! extract_hjs
c-----------------------------------------------------------------------
