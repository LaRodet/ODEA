c***********************************************************************
c                          IO_SPIN_WRITE.F
c***********************************************************************
c Write the spin with respect to time
c
c Input:
c    t         ==> Current time
c    nbod      ==> Number of massive bodies
c    stidal    ==> spin
c    iu        ==> Unit to write to for massive bodies
c    fopenstat ==> The status flag for the open statements
c                    of the output files.

      subroutine io_spin_write(t, nbod, stidal, iu, diro, fopenstat)

      include '../swift.inc'
      include 'io.inc'

c...  Inputs:
      integer iu, nbod
      real*8 t, stidal(NPLMAX)
      character*(*) fopenstat, diro

c...  Internals
      integer ierr, i1st
      data i1st/0/
      save i1st

c----
c...  Executable code

      if(i1st.eq.0) then

         call io_open(iu,trim(diro)//'/'//'spin.out',
     &        fopenstat,'FORMATTED',ierr)
         if(ierr.ne.0) then
            write(*,*) ' ERROR in io_spin_write '
            write(*,*) '  Could not open spin.out '
            call util_exit(1)
         endif

      else

         call io_open(iu,trim(diro)//'/'//'spin.out',
     &        'append','FORMATTED',ierr)

      endif

      write(iu,*) t,stidal(1:nbod)

      close(iu)

      if(i1st.eq.0) then
         i1st=1
      endif

      return
      end   ! io_spin_write
c-----------------------------------------------------------------------
