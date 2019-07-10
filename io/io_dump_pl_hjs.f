c***********************************************************************
c                         IO_DUMP_PL_HJS.F
c***********************************************************************
c Dumps the data for the massive bodies, HJS case
c
c Input:
c    dplfile       ==>  Name of file to write to
c    nbod          ==>  Number of massive bodies
c    oloc          ==>  Link matrix between bodies & orbits
c                      oloc(j,i)=1  : body #i is a satellite in orbit #j
c                      oloc(j,i)=-1 : body #i is a center in orbit #j
c    mass          ==>  Masses of bodies
c    umat          ==>  Conversion matrix Gen. Jacobi => Bary
c    xj, yj, zj    ==>  Initial position in Jacobi coord
c    vxj, vyj, vzj ==>  Initial position in Jacobi coord
c    lclose        ==> .true. --> discard particle if it gets
c                                    too close to a planet.
c    iflgchk       ==>  bit 5 set ==>  include J2 and J4 terms
c    rplsq         ==>  min distance^2 that a tp can get from pl

      subroutine io_dump_pl_hjs(dplfile, nbod, oloc, mass, umat,
     &   xj, yj, zj, vxj, vyj, vzj, lclose, iflgchk, rplsq)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs
      integer nbod, iflgchk, oloc(NPLMAX,NPLMAX)
      real*8 mass(nbod), rplsq(nbod)
      real*8 xj(nbod), yj(nbod), zj(nbod)
      real*8 vxj(nbod), vyj(nbod), vzj(nbod)
      real*8 umat(NPLMAX,NPLMAX)
      character*(*) dplfile
      logical*2 lclose

c...  Internals
      integer j,ierr,i
      real*8 xb(NPLMAX),yb(NPLMAX),zb(NPLMAX)
      real*8 vxb(NPLMAX),vyb(NPLMAX),vzb(NPLMAX)
      real*8 rpl

 123  format(3(1p1e23.16,1x))

c-----
c...  Executable code

      call io_open(7,dplfile,'unknown','formatted',ierr)

      write(7,*) nbod

      call coord_g2b(nbod, umat, mass, xj, yj, zj, vxj, vyj, vzj,
     &   xb, yb, zb, vxb, vyb, vzb)

      do j=1,nbod
         if(lclose) then
            rpl = sqrt(rplsq(j))
            write(7,123) mass(j),rpl
         else
            write(7,123) mass(j)
         endif
         write(7,123) xb(j), yb(j), zb(j)
         write(7,123) vxb(j), vyb(j), vzb(j)
      enddo

      do j = 2,nbod
         write(7,*)(oloc(j,i),i=1,nbod)
      end do

      close(unit = 7)
      return
      end    ! io_dump_pl_hjs.f
c-----------------------------------------------------------------------
