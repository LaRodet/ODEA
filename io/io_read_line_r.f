c***********************************************************************
c                            IO_READ_LINE_R.F
c***********************************************************************
c Read one line from real*8 binary file.
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

      integer function io_read_line_r(iu, id, a, e, inc, capom, omega,
     &   capm)

c...  Inputs:
      integer iu

c...  Output:
      integer id
      real*8 a, e, inc, capom, omega, capm

c...  Internals
      integer ierr

c----
c...  Executable code

      read(iu,iostat=ierr) id, a, e, inc, capom, omega, capm
      io_read_line_r = ierr
      if(ierr.ne.0) then
         return
      endif

      return
      end      ! io_read_line_r
c-----------------------------------------------------------------------
