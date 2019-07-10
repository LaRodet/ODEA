c*************************************************************************
c                            IO_READ_HDR
c*************************************************************************
c Read in header part of the real*4 binary file
c
c Input:
c    iu            ==> unit number to write to
c
c Output:
c    time          ==>  current time (real scalar)
c    nbod          ==>  number of massive bodies (int scalar)
c    nleft         ==>  number of active tp (int scalar)

      integer function io_read_hdr(iu, time, nbod, nleft)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs:
      integer iu

c...  Output
      integer nbod, nleft
      real*8 time

c...  Internals
      real*4 ttmp
      integer*2 nleft2, nbod2
      integer ierr

c----
c...  Executable code


      read(iu,iostat=ierr) ttmp, nbod2, nleft2
      io_read_hdr = ierr
      if(ierr.ne.0) then
         return
      endif

      nbod = nbod2
      nleft = nleft2
      time = ttmp

      return
      end     ! io_read_hdr
c---------------------------------------------------------------------------
