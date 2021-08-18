c***********************************************************************
c                         IO_DUMP_TIDAL.F
c***********************************************************************
c Dumps the data for the massive bodies, HJS case
c
c Input:
c    tidalfile     ==>  Name of file to write to
c    nbod          ==>  Number of massive bodies
c    mass          ==>  Masses of bodies
c    atidal        ==>  Dimensionless moment of inertia
c    qtidal        ==>  Modified tidal factor
c    rtidal        ==>  Radius
c    sx, sy, sz    ==>  Spin components

      subroutine io_dump_tidal(tidalfile, nbod, mass, atidal, qtidal,
     &   rtidal, sx, sy, sz)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs
      integer nbod
      real*8 mass(nbod), atidal(nbod), qtidal(nbod)
      real*8 rtidal(nbod), sx(nbod), sy(nbod), sz(nbod)
      character*(*) tidalfile

c...  Internals
      integer j,ierr

 123  format(3(1p1e23.16,1x))

c-----
c...  Executable code

      call io_open(7,tidalfile,'unknown','formatted',ierr)

      write(7,*) nbod

      do j=1,nbod
         write(7,123) mass(j), atidal(j), qtidal(j), rtidal(j),
     &      sx(j), sy(j), sz(j)
      enddo

      close(unit = 7)
      return
      end    ! io_dump_tidal.f
c-----------------------------------------------------------------------
