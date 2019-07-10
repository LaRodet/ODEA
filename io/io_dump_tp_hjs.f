c***********************************************************************
c                       IO_DUMP_TP_HJS.F
c***********************************************************************
c Dump test particle data, in the HJS case
c
c Inputs:
c    dtpfile          ==>  Name of file to write to
c    nbod             ==>  Number of massive bodies
c    ntp              ==>  number of massive bodies
c    matp             ==> Conversion vectors bary=>Jacobi for tp's
c    xjt, yjt, zjt    ==>  initial position in Jacobi coord
c    vxjt, vyjt, vzjt ==>  initial position in Jacobi coord
c    istat            ==>  status of the test paricles
c                                      istat(i,1) = 0  active
c                                      istat(i,1) = 1 not
c    rstat            ==>  status of the test paricles

      subroutine io_dump_tp_hjs(dtpfile, nbod, ntp, matp,
     &	xjt, yjt, zjt, vxjt, vyjt, vzjt, istat, rstat)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs
      integer ntp, nbod
      real*8 xjt(ntp), yjt(ntp), zjt(ntp)
      real*8 vxjt(ntp), vyjt(ntp), vzjt(ntp)
      real*8 matp(NPLMAX,NTPMAX)
      real*8 rstat(NTPMAX,NSTATR)
      integer istat(NTPMAX,NSTAT)
      integer iflgchk
      character*(*) dtpfile

c...  Internal
      integer i, j, ierr, orbct(NPLMAX)

 123  format(4(1p1e23.16,1x))
 124  format(1(1p1e23.16,1x))

c-----
c...  Executable code

      call io_open(7,dtpfile,'unknown','formatted',ierr)

      write(7,*) ntp
      do i=1,ntp
         do j=1,nbod
            orbct(j) = 0
         end do
         do j = 1,nbod
            if (matp(j,i).ne.0.0d0) orbct(j) = -1
         end do
         write(7,123) xjt(i),yjt(i),zjt(i)
         write(7,123) vxjt(i),vyjt(i),vzjt(i)
         write(7,*) (orbct(j),j=1,nbod)
         write(7,*) (istat(i,j),j=1,NSTAT)
         write(7,123) (rstat(i,j),j=1,NSTATR)
      end do

      close(unit = 7)

      return
      end    ! io_dump_tp_hjs
c-----------------------------------------------------------------------
