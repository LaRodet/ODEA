c***********************************************************************
c                            IO_WRITE_HDR_ODEA
c***********************************************************************
c Write out header part of the binary file
c
c Input:
c    iu    ==> unit number to write to
c    time  ==> current time
c    nbod  ==> number of massive bodies
c    ntp   ==> number of massive bodies
c    istat ==> status of the test paricles

      subroutine io_write_hdr_odea(iu, time, nbod, ntp, istat)

      include '../swift.inc'

c...  Inputs:
      integer nbod, ntp, istat(NTPMAX,NSTAT), iu
      real*8 time

c...  Internals
      integer i, nleft

c----
c...  Executable code


c...  calculate number of remaining test particles
      nleft = 0
      do i=1,ntp
         if(istat(i,1).eq.0) then
            nleft = nleft + 1
         endif
      enddo

      write(iu) time, nbod, nleft

      return
      end     ! io_write_hdr_odea.f
c-----------------------------------------------------------------------
