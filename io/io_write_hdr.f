c***********************************************************************
c                            IO_WRITE_HDR
c***********************************************************************
c Write out header part of the real*4 binary file
c
c Input:
c    iu    ==> unit number to write to
c    time  ==> current time
c    nbod  ==> number of massive bodies
c    ntp   ==> number of massive bodies
c    istat ==> status of the test particles

      subroutine io_write_hdr(iu,time,nbod,ntp,istat)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs:
      integer nbod, ntp, istat(NTPMAX,NSTAT), iu
      real*8 time

c...  Internals
      integer i
      real*4 ttmp
      integer*2 nleft, nbod2

c----
c...  Executable code

c...  Calculate number of remaining test particles
      nleft = 0
      do i=1,ntp
         if(istat(i,1).eq.0) then
            nleft = nleft + 1
         endif
      enddo

      nbod2 = nbod

      ttmp = time

      write(iu) ttmp, nbod2, nleft

      return
      end     ! io_write_hdr
c-----------------------------------------------------------------------
