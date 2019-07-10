c***********************************************************************
c                          IO_ENERGY_WRITE_ODEA.F
c***********************************************************************
c Does the write for anal_jacobi_write_odea
c
c Input:
c    i1st           ==>  =0 if first write, =1 if not
c    t              ==>  Current time
c    energy         ==>  Total energy
c    eltot          ==>  Components of total angular momentum
c    dt             ==>  Time step
c    iu             ==>  Unit to write to
c    fopenstat      ==>  The status flag for the open statements
c                         of the output files.

      subroutine io_energy_write_odea(i1st, t, energy, eltot, dt, iu,
     &   diro, fopenstat)

      include '../swift.inc'

c...  Inputs:
      integer iu, i1st
      real*8 t, energy, eltot(3), dt
      character*(*) fopenstat, diro

c...  Internals
      integer ierr

 2    format(1x, 1p1e12.5, 4(2x,1p1e23.16), 1x, 1p1e12.5)

c----
c...  Executable code

      if(i1st.eq.0) then

         call io_open(iu, trim(diro)//'/'//'energy.out',
     &      fopenstat, 'FORMATTED', ierr)
         if(ierr.ne.0) then
            write(*,*) ' Error in anal_energy_write_odea: '
            write(*,*) ' Could not open energy.out '
            call util_exit(1)
         endif

      else

         call io_open(iu, trim(diro)//'/'//'energy.out', 'append',
     &      'FORMATTED', ierr)

      endif

      write(iu,2) t, energy, eltot, dt

      close(iu)

      return
      end ! io_energy_write_odea
c-----------------------------------------------------------------------
