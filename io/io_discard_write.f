c***********************************************************************
c                            IO_DISCARD_WRITE
c***********************************************************************
c Checks to see if a particle was removed since te last time
c  the subroutine was called.
c If so write out the position and velocity of tp and planets.
c
c Input:
c    init          ==>  Initialize flag if = 0 initialize and return
c                                          = 1 run through
c    time          ==>  Current time
c    nbod          ==>  Number of massive bodies
c    ntp           ==>  Number of test bodies
c    x, y, z       ==>  Current position (Jacobi or heliocentric)
c    vx, vy, vz    ==>  Current velocity (Jacobi or heliocentric)
c    xt, yt, zt    ==>  Current tp position (Jacobi or heliocentric)
c    vxt, vyt, vzt ==>  Current tp velocity (Jacobi or heliocentric)
c    istat         ==>  Status of the test paricles
c                         istat(i,1) = 0 ==> active:  = 1 not
c                         istat(i,2) = -1 ==> Danby did not work
c    rstat         ==>  status of the test paricles
c                         rstat(i,1) time of discard.
c                         rstat(i,2) closest approach to a planet
c    iu            ==> unit number to write to
c    rname         ==> output file name
c    fopenstat     ==>  The status flag for the open statements
c                          of the output files.
c    nleft         ==>  Number of active test bodies

      subroutine io_discard_write(init, time, nbod, ntp, x, y, z,
     &   vx, vy, vz, xt, yt, zt, vxt, vyt, vzt, istat,
     &   rstat, iu, rname, fopenstat, nleft)

      include '../swift.inc'

c...  Inputs:
      integer init, nbod, ntp, iu
      real*8 time
      real*8 rstat(NTPMAX,NSTATR)
      integer istat(NTPMAX,NSTAT)
      real*8 x(nbod), y(nbod), z(nbod)
      real*8 vx(nbod), vy(nbod), vz(nbod)
      real*8 xt(ntp), yt(ntp), zt(ntp)
      real*8 vxt(ntp), vyt(ntp), vzt(ntp)
      character*(*) rname, fopenstat

c...  Outputs:
      integer nleft

c...  Internals
      integer istold(NTPMAX), i, in, irem, j
      integer iwrite, ierr, izero(NSTAT)
      real*8 rzero(NSTATR)

      save istold, izero, iwrite, rzero

 1000 format(1x, 1p1e23.16)
 2000 format(3x, i5, 100(1x,i5))
 3000 format(4(1p1e23.16,1x))
c----
c...  Executable code

c...  If initialize flag=0
      if(init.eq.0) then

         do i=1,ntp
            istold(i) = istat(i,1)
         enddo

         do i=1,NSTAT
            izero(i) = 0
         enddo
         do i=1,NSTATR
            rzero(i) = 0.0d0
         enddo
         iwrite = 0
         return
      endif

c...  If initialize flag=1
      irem = 0
      nleft = 0
      do i=1,ntp
         if(istat(i,1).eq.0) then
            nleft = nleft + 1
         endif
         if(istold(i).ne.istat(i,1)) then
            if(irem.eq.0) then
               if(iwrite.eq.0) then
                  call io_open(iu, rname, fopenstat, 'FORMATTED', ierr)

c...              If there was an error and fopenstat='append' then
c...              try to open as new
                  if(ierr.ne.0) then
                     if( (fopenstat(1:6).eq.'append') .or.
     &                    (fopenstat(1:6).eq.'APPEND') ) then
                        call io_open(iu, rname, 'new', 'FORMATTED',
     &                     ierr)
                     endif
                  endif

                  if(ierr.ne.0) then
                     write(*,*) ' ERROR in io_discard_write: '
                     write(*,*) '  Could not open discard output file'
                     call util_exit(1)
                  endif
               else
                  call io_open(iu, rname, 'append', 'FORMATTED', ierr)
               endif
               write(iu,1000) time
            endif
            iwrite = 1
            irem = irem + 1
            write(iu,2000) i, (istat(i,j),j=1,NSTAT)
            write(iu,3000) (rstat(i,j),j=1,NSTATR)
            write(iu,3000) xt(i), yt(i), zt(i)
            write(iu,3000) vxt(i), vyt(i), vzt(i)
         endif
         istold(i) = istat(i,1)
      enddo

      if(irem.ne.0) then
         do i=2,nbod
            in = -1*i
            write(iu,2000) in, (izero(j),j=1,NSTAT)
            write(iu,3000) (rzero(j),j=1,NSTATR)
            write(iu,3000) x(i), y(i), z(i)
            write(iu,3000) vx(i), vy(i), vz(i)
         enddo
         in = 0
         write(iu,2000) in, (izero(j),j=1,NSTAT)
      endif

      close(iu)

      return
      end     !io_discard_write
c-----------------------------------------------------------------------
