c***********************************************************************
c                            IO_WRITE_FRAME_ODEA
c***********************************************************************
c Write out a whole frame to a binary file for
c both massive and test particles (Generalized Jacobi coordinates)
c
c Inputs:
c    time             ==> Current time
c    nbod             ==> Number of massive bodies
c    ntp              ==> Number of test particles
c    oloc             ==> Link between bodies & orbits
c    mass             ==> Mass of bodies
c    eta, mu          ==> Masses of centers & satellites
c    umat             ==> Conversion matrix Jacobi => Barycentric
c    xj, yj, zj       ==> Current position in Jacobi coord
c    vxj, vyj, vzj    ==> Current velocity in Jacobi coord
c    matp             ==> Conversion vectors for tp's
c    etatp            ==> Masses of centers for tp's
c    xjt, yjt, zjt    ==> Current tp position in Jacobi coord
c    vxjt, vyjt, vzjt ==> Current tp velocity in Jacobi coord
c    istat            ==> Status of the test paricles
c    oname            ==> Output file name
c    iu               ==> Unit number to write to
c    fopenstat        ==> Status flag for the open
c                          statements of the output files.
c    iflgchk          ==> First two bits indicate real*4 or *8

      subroutine io_write_frame_odea(time, nbod, ntp, oloc, matp, umatp,
     &   mat, umat, mass, eta, mu, xj, yj, zj, vxj, vyj, vzj, etatp,
     &   xjt, yjt, zjt, vxjt, vyjt, vzjt, istat, oname, iu, fopenstat,
     &   iflgchk, matr)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs:
      integer nbod, ntp, iu
      real*8 mass(nbod), eta(nbod), mu(nbod),time
      integer istat(ntp), oloc(NPLMAX,NPLMAX), iflgchk
      real*8 matp(NPLMAX,NTPMAX), umat(NPLMAX,NPLMAX), etatp(NTPMAX)
      real*8 mat(NPLMAX,NPLMAX), umatp(NPLMAX,NTPMAX)
      real*8 xj(nbod), yj(nbod), zj(nbod)
      real*8 vxj(nbod), vyj(nbod), vzj(nbod)
      real*8 xjt(ntp), yjt(ntp), zjt(ntp)
      real*8 vxjt(ntp), vyjt(ntp), vzjt(ntp)
      real*8 matr(3,3)
      character*(*) oname, fopenstat

c...  Internals
      integer i, id, j, k
      integer ialpha, ierr
      logical ok
      real*8 mat0(NPLMAX,NPLMAX), matp0(NPLMAX,NTPMAX), etatp0(NTPMAX)
      real*8 xj0(nbod), yj0(nbod), zj0(nbod)
      real*8 vxj0(nbod), vyj0(nbod), vzj0(nbod)
      real*8 xjt0, yjt0, zjt0, vxjt0, vyjt0, vzjt0
      real*8 xt, yt, zt, vxt, vyt, vzt
      real*8 xb(nbod), yb(nbod), zb(nbod)
      real*8 vxb(nbod), vyb(nbod), vzb(nbod)
      real*8 xbt, ybt, zbt, vxbt, vybt, vzbt
      real*8 a, e, inc, capom, omega, capm
      real*8 gm, eta0(NPLMAX), mu0(NPLMAX)
      integer i1st    ! =0 first time through; =1 after
      data i1st/0/
      save i1st, mat0, eta0, mu0, matp0, etatp0

c----
c...  Executable code

c...  If first time through open file
      if(i1st.eq.0) then

         call io_open(iu,oname,fopenstat,'UNFORMATTED',ierr)
         if(ierr.ne.0) then
            write(*,*) ' SWIFT ERROR: in io_write_frame: '
            write(*,*) '     Could not open binary output file:'
            call util_exit(1)
         endif

         mat0 = mat
         mu0 = mu
         eta0 = eta

         matp0 = matp
         etatp0 = etatp

         i1st = 1
      else
         call io_open(iu, oname, 'append', 'UNFORMATTED', ierr)
      endif

      if (btest(iflgchk,0)) then
         call io_write_hdr(iu, time, nbod, ntp, istat)
      else
         call io_write_hdr_r(iu, time, nbod, ntp, istat)
      end if

      call coord_g2b(nbod, umat, mass, xj, yj, zj, vxj, vyj, vzj,
     &   xb, yb, zb, vxb, vyb, vzb)

      call coord_b2g(nbod, mat0, mass, xb, yb, zb, vxb, vyb, vzb,
     &   xj0, yj0, zj0, vxj0, vyj0, vzj0)

c...  Write out massive bodies
      do i=2,nbod
         gm = eta0(i) + mu0(i)
         id = -1*i
         call orbel_xv2el(xj0(i), yj0(i), zj0(i), vxj0(i), vyj0(i),
     &      vzj0(i), gm, ialpha, a, e, inc, capom, omega, capm)
         if (btest(iflgchk,0)) then
            call io_write_line(iu, id, a, e, inc, capom, omega, capm)
         else
            call io_write_line_r(iu, id, a, e, inc, capom, omega, capm)
         end if
      end do

c...  write out test particles.
      do i=1,ntp
         if (istat(i).eq.0) then

            call coord_g2b_tp(nbod, umatp(:,i), xj, yj, zj, vxj, vyj,
     &         vzj, xjt(i), yjt(i), zjt(i), vxjt(i), vyjt(i), vzjt(i),
     &         xbt, ybt, zbt, vxbt, vybt, vzbt)

            call coord_vb2vg_tp(nbod, matp0(:,i), xb, yb, zb,
     &         xbt, ybt, zbt, xjt0, yjt0, zjt0)

            call coord_vb2vg_tp(nbod, matp0(:,i), vxb, vyb, vzb,
     &         vxbt, vybt, vzbt, vxjt0, vyjt0, vzjt0)

            xt = matr(1,1)*xjt0+matr(1,2)*yjt0+matr(1,3)*zjt0
            yt = matr(2,1)*xjt0+matr(2,2)*yjt0+matr(2,3)*zjt0
            zt = matr(3,1)*xjt0+matr(3,2)*yjt0+matr(3,3)*zjt0
            vxt = matr(1,1)*vxjt0+matr(1,2)*vyjt0+matr(1,3)*vzjt0
            vyt = matr(2,1)*vxjt0+matr(2,2)*vyjt0+matr(2,3)*vzjt0
            vzt = matr(3,1)*vxjt0+matr(3,2)*vyjt0+matr(3,3)*vzjt0

            call orbel_xv2el(xt, yt, zt, vxt, vyt, vzt,
     &         etatp0(i), ialpha, a, e, inc, capom, omega, capm)

            if (btest(iflgchk,0)) then
               call io_write_line(iu, i, a, e, inc, capom, omega, capm)
            else
               call io_write_line_r(iu, i, a, e, inc, capom, omega,capm)
            end if

         endif
      enddo

      close(iu)
      return
      end      ! io_write_frame_odea
c-----------------------------------------------------------------------
