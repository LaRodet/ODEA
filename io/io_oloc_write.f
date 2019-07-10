c***********************************************************************
c                          IO_OLOC_WRITE.F
c***********************************************************************
c Write the hierarchy with respect from time
c
c Input:
c    t         ==> Current time
c    oloc      ==> Link matrix between bodies & orbits
c                   oloc(j,i)=1 : body #i is a satellite in orbit #j
c                   oloc(j,i)=-1: body #i is a center in orbit #j
c    oloctp    ==> Link matrix between bodies & tp orbits
c                   oloc(j,i)=1: body #j is a satellite in tp #i's orbit
c                   oloc(j,i)=-1: body #j is a center in tp #i's orbit
c    iu        ==> Unit to write to for massive bodies
c    iutp      ==> Unit to write to for tp
c    fopenstat ==> The status flag for the open statements
c                    of the output files.

      subroutine io_oloc_write(t, nbod, ntp, oloc, oloct, iu, iutp,
     &   diro, fopenstat)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs:
      integer iu, iutp, nbod, oloc(NPLMAX,NPLMAX)
      integer ntp, oloct(NPLMAX,NTPMAX)
      real*8 t
      character*(*) fopenstat, diro

c...  Internals
      integer ierr, i1st
      data i1st/0/
      save i1st

c----
c...  Executable code

      if(i1st.eq.0) then

         call io_open(iu,trim(diro)//'/'//'oloc.out',
     &        fopenstat,'FORMATTED',ierr)
         if(ierr.ne.0) then
            write(*,*) ' ERROR in anal_oloc_write '
            write(*,*) '  Could not open oloc.out '
            call util_exit(1)
         endif

      else

         call io_open(iu,trim(diro)//'/'//'oloc.out',
     &        'append','FORMATTED',ierr)

      endif

      write(iu,*) t,oloc(2:nbod,1:nbod)

      close(iu)

      if(i1st.eq.0) then

         call io_open(iutp, trim(diro)//'/'//'oloc_tp.out',
     &        fopenstat, 'FORMATTED', ierr)
         if(ierr.ne.0) then
            write(*,*) ' ERROR in io_oloc_write '
            write(*,*) '  Could not open oloc_tp.out '
            call util_exit(1)
         endif

      else

         call io_open(iutp, trim(diro)//'/'//'oloc_tp.out',
     &        'append', 'FORMATTED', ierr)

      endif

      write(iutp,*) t, oloct(2:nbod,1:ntp)

      close(iutp)

      if(i1st.eq.0) then
         i1st=1
      endif

      return
      end                       ! io_oloc_write
c-----------------------------------------------------------------------
