c***********************************************************************
c                            IO_WRITE_LINE_R
c***********************************************************************
c Write out one line to real*8 binary file.
c
c Inputs:
c    iu     ==> unit number to write to
c    a      ==> semi-major axis or pericentric distance if a parabola
c    e      ==> eccentricity
c    inc    ==> inclination
c    capom  ==> longitude of ascending node
c    omega  ==> argument of perihelion
c    capm   ==> mean anomaly

      subroutine io_write_line_r(iu, id, a, e, inc, capom, omega, capm)

      include '../swift.inc'

c...  Inputs:
      integer iu, id
      real*8 a, e, inc, capom, omega, capm

c----
c...  Executable code

      write(iu) id, a, e, inc, capom, omega, capm

      return
      end      ! io_write_line_r
c-----------------------------------------------------------------------
