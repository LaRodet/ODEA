c***********************************************************************
c                              IO_INIT_TIDAL.F
c***********************************************************************
c Reads in the tidal data for the massive bodies
c
c Input:
c    tidalfile     ==>  Name of file to write to
c
c Output:
c    nbod          ==>  Number of massive bodies
c    mass          ==>  Masses of bodies
c    atidal        ==>  Dimensionless moment of inertia
c    qtidal        ==>  Modified tidal factor
c    rtidal        ==>  Radius
c    stidal        ==>  Spin

      subroutine io_init_tidal(tidalfile, mass, atidal, qtidal, rtidal,
     &   stidal)

      include '../swift.inc'
      include 'io.inc'

c...  Input
      character*(*) tidalfile

c...  Output
      integer nbod
      real*8 mass(NPLMAX), atidal(NPLMAX), qtidal(NPLMAX)
      real*8 rtidal(NPLMAX), stidal(NPLMAX)

c...  Internal
      integer j,ierr
      real*8 rpl, mtot, vsat, vcen
      real*8 xb(NPLMAX), yb(NPLMAX), zb(NPLMAX)
      real*8 vxb(NPLMAX), vyb(NPLMAX), vzb(NPLMAX)

c-----
c...  Executable code

      write(*,*) 'Tidal data file is ', tidalfile
      call io_open(7, tidalfile, 'old', 'formatted', ierr)

c...  Read number of planets
      read(7,*) nbod

      if(nbod.gt.NPLMAX) then
         write(*,*) ' SWIFT ERROR: in io_init_pl: '
         write(*,*) '   The number of massive bodies,',nbod,','
         write(*,*) '   is too large, it must be less than',NPLMAX
         call util_exit(1)
      endif

      write(*,23) nbod
 23   format(/,'Number of bodies is ',i3,/,
     &   "For each, list mass, alpha, Q', R, spin")

c... Read information relative to bodies: masses, bary pos & vels.
      do j=1,nbod
          read(7,*) mass(j), atidal(j), qtidal(j), rtidal(j), stidal(j)
          write(*,*) mass(j), atidal(j), qtidal(j), rtidal(j), stidal(j)
      enddo

      return

      end     ! io_init_tidal.f
c-----------------------------------------------------------------------
