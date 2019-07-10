c***********************************************************************
c                            IO_WRITE_LINE
c***********************************************************************
c Write out one line to real*4 binary file.
c
c Input:
c    iu    ==> unit number to write to
c    a     ==> semi-major axis or pericentric distance if a parabola
c    e     ==> eccentricity
c    inc   ==> inclination
c    capom ==> longitude of ascending node
c    omega ==> argument of perihelion
c    capm  ==> mean anomaly

      subroutine io_write_line(iu, id, a, e, inc, capom, omega, capm)

      include '../swift.inc'

c...  Inputs:
      integer iu,id
      real*8 a, e, inc, capom, omega, capm

c...  Internals
      integer*2 id2
      real*4 a4, e4, inc4, capom4, omega4, capm4

c----
c...  Executable code

      id2 = id

      a4 = a
      e4 = e
      inc4 = inc
      capom4 = capom
      capm4 = capm
      omega4 = omega

      write(iu) id2, a4, e4, inc4, capom4, omega4, capm4

      return
      end      ! io_write_line_r
c-----------------------------------------------------------------------
