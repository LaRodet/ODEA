c***********************************************************************
c                            UTIL_EXIT.F
c***********************************************************************
c Exits program
c
c Input:
c    iflg ==> status of exit
c                = 0 if normal exit
c                = 1 if exit because error

      subroutine util_exit(iflg)

      include '../swift.inc'

c...  Inputs:
      integer iflg

c-----
c...  Executable code

      write(*,*) ' '

      if(iflg.eq.0) then
         write(*,*) 'Normal termination'
      else
         write(*,*) 'Terminating due to ERROR!'
      endif

      write(*,*) '----------------------------------------------------'

      stop
      end  ! util_exit
c-----------------------------------------------------------------------
