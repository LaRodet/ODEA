c***********************************************************************
c                            IO_OPEN.F
c***********************************************************************
c Open files
c
c Input:
c    iu          ==>  unit number
c    fname       ==>  file name
c    fopenstat   ==>  The status flag for the open
c                      statements of the output files.
c    format      ==>  format string
c
c Output:
c    ierr        ==>  output from iostat

      subroutine io_open(iu, fname, fopenstat, format, ierr)

      include '../swift.inc'

c...  Inputs:
      integer iu
      character*(*) fname,fopenstat,format

c...  Outputs:
      integer ierr

c----
c...  Executable code

      if( (fopenstat(1:6).eq.'append') .or.
     &   (fopenstat(1:6).eq.'APPEND') ) then
         open(unit=iu, file=fname, status='old', access='append',
     &      form=format,iostat=ierr)
         if(ierr.ne.0) then
            write(*,*) 'Warning:  Could not open ',fname,' with'
            write(*,*) '          position=append.'
            write(*,*) '          access=append.'
            write(*,*) '          Will open as status=new'
            open(unit=iu, file=fname, status='new',
     &         form=format, iostat=ierr)
         endif
      else
         open(unit=iu, file=fname, status=fopenstat,
     &        form=format, iostat=ierr)

      endif

      return
      end   ! io_open
c-----------------------------------------------------------------------
