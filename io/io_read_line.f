c***********************************************************************
c                            IO_READ_LINE.F
c***********************************************************************
c Read one line from real*4 binary file.
c
c Input:
c    iu       ==> Unit number to read to
c
c Outputs:
c    a        ==> Semi-major axis or pericentric distance if a parabola
c    e        ==> Eccentricity
c    inc      ==> Inclination
c    capom    ==> Longitude of ascending node
c    omega    ==> Argument of perihelion
c    capm     ==> Mean anomaly
c
c Returns:
c    io_read_line_r    ==>   =0 read ok else failed

      integer function io_read_line(iu, id, a, e, inc, capom, omega,
     &   capm)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs:
      integer iu

c...  Outputs:
      integer id
      real*8 a, e, inc, capom, omega, capm

c...  Internals
      integer*2 id2
      real*4 a4, e4, inc4, capom4, omega4, capm4
      integer ierr

c----
c...  Executable code

      read(iu,iostat=ierr) id2, a4, e4, inc4, capom4, omega4, capm4
      io_read_line = ierr
      if(ierr.ne.0) then
         return
      endif

      id = id2

      a = a4
      e = e4
      inc = inc4
      capom = capom4
      capm = capm4
      omega = omega4

      return
      end      ! io_read_line
c-----------------------------------------------------------------------
