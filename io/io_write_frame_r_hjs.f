c***********************************************************************
c                            IO_WRITE_FRAME_R_HJS
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

      subroutine io_write_frame_r_hjs(time, nbod, ntp, oloc, matp, umat,
     &   mass, eta, mu, xj, yj, zj, vxj, vyj, vzj, etatp, xjt, yjt, zjt,
     &   vxjt, vyjt, vzjt, istat, oname, iu, fopenstat, matr)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs:
      integer nbod, ntp, iu
      real*8 mass(nbod), eta(nbod), mu(nbod), time
      integer istat(ntp), oloc(NPLMAX,NPLMAX)
      real*8 matp(NPLMAX,NTPMAX), umat(NPLMAX,NPLMAX), etatp(ntp)
      real*8 xj(nbod), yj(nbod), zj(nbod)
      real*8 vxj(nbod), vyj(nbod), vzj(nbod)
      real*8 xjt(ntp), yjt(ntp), zjt(ntp)
      real*8 vxjt(ntp), vyjt(ntp), vzjt(ntp)
      real*8 matr(3,3)
      character*(*) oname, fopenstat

c...  Internals
      integer i, id
      integer ialpha, ierr
      logical ok
      real*8 xt, yt, zt, vxt, vyt, vzt
      real*8 gm, a, e, inc, capom, omega, capm
      integer i1st    ! =0 first time through; =1 after
      data i1st/0/
      save i1st

c----
c...  Executable code

c...  if first time through open file
      if(i1st.eq.0) then
         call io_open(iu,oname,fopenstat,'UNFORMATTED',ierr)
         if(ierr.ne.0) then
            write(*,*) ' SWIFT ERROR: in io_write_frame: '
            write(*,*) '     Could not open binary output file:'
            call util_exit(1)
         endif
         i1st = 1
      else
         call io_open(iu, oname, 'append', 'UNFORMATTED', ierr)
      endif

      call io_write_hdr_r(iu, time, nbod, ntp, istat)

c...  write out planets
      do i=2,nbod
         gm = eta(i) + mu(i)
         id = -1*i
         call orbel_xv2el(xj(i), yj(i), zj(i), vxj(i), vyj(i), vzj(i),
     &      gm, ialpha, a, e, inc, capom, omega, capm)
         call io_write_line_r(iu, id, a, e, inc, capom, omega, capm)
      end do

c...  write out test particles.
      do i=1,ntp
         if (istat(i).eq.0) then

            xt = matr(1,1)*xjt(i)+matr(1,2)*yjt(i)+matr(1,3)*zjt(i)
            yt = matr(2,1)*xjt(i)+matr(2,2)*yjt(i)+matr(2,3)*zjt(i)
            zt = matr(3,1)*xjt(i)+matr(3,2)*yjt(i)+matr(3,3)*zjt(i)
            vxt = matr(1,1)*vxjt(i)+matr(1,2)*vyjt(i)+matr(1,3)*vzjt(i)
            vyt = matr(2,1)*vxjt(i)+matr(2,2)*vyjt(i)+matr(2,3)*vzjt(i)
            vzt = matr(3,1)*vxjt(i)+matr(3,2)*vyjt(i)+matr(3,3)*vzjt(i)

            call orbel_xv2el(xt, yt, zt, vxt, vyt, vzt,
     &         etatp(i), ialpha, a, e, inc, capom, omega, capm)
            call io_write_line_r(iu, i, a, e, inc, capom, omega, capm)
         endif
      enddo

      close(iu)
      return
      end      ! io_write_frame_r_hjs
c-----------------------------------------------------------------------
