c***********************************************************************
c                            IO_READ_HDR_R.F
c***********************************************************************
c Read in header part of the real*8 binary file
c
c Input:
c    iu            ==> unit number to write to
c
c Output:
c    time          ==>  current time (real scalar)
c    nbod          ==>  number of massive bodies (int scalar)
c    nleft         ==>  number of active tp (int scalar)
c
c Returns:
c    io_read_hdr_r   ==>   =0 read ok, fails otherwise

      integer function io_read_hdr_r(iu, time, nbod, nleft)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs:
      integer iu

c...  Output
      integer nbod, nleft
      real*8 time

c...  Internals
      integer ierr

c----
c...  Executable code


      read(iu,iostat=ierr) time, nbod, nleft
      io_read_hdr_r = ierr
      if(ierr.ne.0) then
         return
      endif

      return
      end     ! io_read_hdr_r
c-----------------------------------------------------------------------
