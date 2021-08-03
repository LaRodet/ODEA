c***********************************************************************
c                             ODEA.F
c***********************************************************************
c  This code is a realization of Hierarchical Jacobi Symplectic
c  N-body mapping method for hierarchical stellar systems in evolving
c  architecture (Rodet et al. 2019).
c  It is based on SWIFT HJS (Beust 2003) and on the SWIFT implementation
c  of Levison and Duncan (1994)
c
c  To run, need 3 input files. The code prompts for the file names,
c  but examples are :
c     parameter file ==> param.in
c     planet file ==> pl.in
c     test particle file ==> tp.in

      include '../swift.inc'

      real*8, parameter :: dr = 1.7453292519943294d-2

      real*8 xjt(NTPMAX), yjt(NTPMAX), zjt(NTPMAX)
      real*8 vxjt(NTPMAX), vyjt(NTPMAX), vzjt(NTPMAX)

      real*8 mass(NPLMAX)
      real*8 xj(NPLMAX), yj(NPLMAX), zj(NPLMAX)
      real*8 vxj(NPLMAX), vyj(NPLMAX), vzj(NPLMAX)

      integer istat(NTPMAX,NSTAT), i1st, istati(NSTAT)
      real*8 rstat(NTPMAX,NSTATR),rstati(NSTATR)
      integer oloc(NPLMAX,NPLMAX), oloct(NPLMAX,NTPMAX)
      integer nbod, ntp, nleft
      integer iflgchk, iuf, iud, iue, iuo, iuotp, i, j, ierr
      real*8 umat(NPLMAX,NPLMAX), mat(NPLMAX,NPLMAX)
      real*8 umatp(NPLMAX,NTPMAX), matp(NPLMAX,NTPMAX)
      real*8 mu(NPLMAX), eta(NPLMAX), etatp(NTPMAX)
      real*8 matr(3,3)

      real*8 t0, tstop, dt, dtout, dtdump
      real*8 t, tout, tdump, tfrac, eoff

      real*8 rmin, rmax, rmaxu, qmin, rplsq(NPLMAX) ! Ignored
      real*8 energy, energy0
      logical*2 lclose

      character*80 outfile, inparfile, inplfile, intpfile, fopenstat
      character*80 diro, dirs, gname, dataname, genfile

      integer ialpha, nsta, fverb, iverb
      real*8 a, e, inc, capom, omega, capm
      logical ok

 998  format(' Time = ',1p1e12.5,': fraction done = ',0pf5.3,
     &           ': Number of active tp =',i4)
 999  format(a)

c-----
c...  Executable code

c...  Beginning the simulation
      call util_version_odea

c..   Get gen file for copy in the simulation directory
      write(*,*) 'Enter name of generic file:'
      read(*,999) genfile

c...  Get data for the run and the test particles
      write(*,*) 'Enter name of parameter data file: '
      read(*,999) inparfile
      call io_init_param(inparfile, t0, tstop, dt, dtout, dtdump,
     &   iflgchk, rmin, rmax, rmaxu, qmin, lclose, diro, dirs,
     &   gname, outfile, fopenstat)

c...  Prompt and read name of massive bodies data file
      write(*,*) ' '
      write(*,*) 'Enter name of massive bodies data file: '
      read(*,999) inplfile
      call io_init_pl_hjs(inplfile, lclose, iflgchk, nbod, oloc,
     &   mass, eta, mu, mat, umat, xj, yj, zj, vxj, vyj, vzj, rplsq)

c...  Get data for the run and the test particles
      write(*,*) 'Enter name of test particle data file: '
      read(*,999) intpfile
      call io_init_tp_hjs(intpfile, nbod, ntp, oloc, oloct, mass, umat,
     &   etatp, matp, umatp, xjt, yjt, zjt, vxjt, vyjt, vzjt, istat,
     &   rstat)

      write(*,*) 'Enter verbose frequency: '
      read(*,*) fverb
      write(*,*) ' fverb = ',fverb

c... Copy initial files input into work directory
      if ((fopenstat(1:3).eq.'new')
     &   .or.(fopenstat(1:3).eq.'NEW')) then
         dataname = 'cp '//trim(inparfile)//' '//trim(diro)
         call system(dataname)
         dataname = 'cp '//trim(inplfile)//' '//trim(diro)
         call system(dataname)
         dataname = 'cp '//trim(intpfile)//' '//trim(diro)
         call system(dataname)
         dataname = 'cp '//trim(genfile)//' '//trim(diro)
         call system(dataname)
         if (ntp.gt.0) then
            dataname = 'cp matpass.dat '//trim(diro)
            call system(dataname)
         end if
      end if

c...  Read the matrix that encodes the disk plane
      if (ntp.gt.0) then

         open(unit=17,file=trim(diro)//'/'//'matpass.dat',status='old',
     &        iostat=ierr)
         if(ierr.ne.0) then
            write(*,*) 'Error:  Could not open matpass.dat'
            call util_version_odea(1)
         endif

         read(17,*)matr(1,1), matr(1,2), matr(1,3)
         read(17,*)matr(2,1), matr(2,2), matr(2,3)
         read(17,*)matr(3,1), matr(3,2), matr(3,3)
         close(17)

      end if

c...  Initialize initial time and times for first output and first dump
      t = t0
      tout = t0 + dtout
      tdump = t0 + dtdump

      iuf = 20
      iud = 40
      iue = 60
      iuo = 80
      iuotp = 100

c...   Do the initial io write
      if(btest(iflgchk,0).or.btest(iflgchk,1))  then ! bit 0 or 1 is set
         call io_write_frame_odea(t0, nbod, ntp, oloc, matp, umatp, mat,
     &      umat, mass, eta, mu, xj, yj, zj, vxj, vyj, vzj, etatp,
     &      xjt, yjt, zjt, vxjt, vyjt, vzjt, istat(:,1),
     &      trim(diro)//'/'//outfile, iuf, fopenstat, iflgchk, matr)
      endif

      if(btest(iflgchk,2))  then ! bit 2 is set
         eoff = 0.0d0
         call anal_energy_write_odea(t, nbod, umat, mass, xj, yj, zj,
     &      vxj, vyj, vzj, iue, fopenstat, diro, eoff, energy0, dt)
      endif

      if(btest(iflgchk,3))  then ! bit 3 is set
         print*, "No Jacobi integral in HJS yet"
      end if

c...  Initialize discard io routine
      if(btest(iflgchk,4))  then ! bit 4 is set
         call io_discard_write(0, t, nbod, ntp, xj, yj, zj, vxj, vyj,
     &      vzj, xjt, yjt, zjt, vxjt, vyjt, vzjt, istat, rstat, iud,
     &      trim(diro)//'/'//'discard.out', fopenstat, nleft)
      endif

      call io_oloc_write(t,nbod,ntp,oloc,oloct,iuo,iuotp,diro,fopenstat)

      nleft = ntp
      i1st = 0
      nsta = 1
      iverb = 0

      write(*,*) ' ************** MAIN LOOP ****************** '

      ok = .true.

      do while (ok)
         call step_kdk_odea(i1st, t, nbod, ntp, oloc, oloct, mat, umat,
     &      matp, umatp, mass, eta, mu, etatp, xj, yj, zj, vxj, vyj,
     &      vzj, xjt, yjt, zjt, vxjt, vyjt, vzjt, istat, rstat, dt,
     &      iflgchk)

         t = t + nsta*dt
c------------

         if (t.ge.tout) then
            iverb = iverb+1
            if (mod(iverb,fverb).eq.0) then

               print*,'t=',t

               do j=2,nbod
                  call orbel_xv2el(xj(j), yj(j), zj(j), vxj(j), vyj(j),
     &               vzj(j), eta(j)+mu(j), ialpha, a, e, inc, capom,
     &               omega, capm)
                  print*, j, sngl(a), sngl(e), sngl(inc)/dr
               end do

            end if
         end if

         ok = (ok.and.(t.le.tstop))

c------------

         if (btest(iflgchk,4))  then ! bit 4 is set
            do i=1,ntp
               if (istat(i,1).eq.0) THEN
                  do j=1,NSTAT
                     istati(j) = istat(i,j)
                     rstati(j) = rstat(i,j)
                  end do
                  call discard_hjs(t, dt, nbod, i, mass, xj, yj, zj,
     &               vxj, vyj, vzj, xjt(i), yjt(i), zjt(i), vxjt(i),
     &               vyjt(i), vzjt(i), umatp(:,i), etatp(i), rmin, rmax,
     &               rmaxu, qmin, lclose, rplsq, istati, rstati)
                  do j=1,NSTAT
                     istat(i,j) = istati(j)
                     rstat(i,j) = rstati(j)
                  end do
               end if
            end do
            call io_discard_write(1, t, nbod, ntp, xj, yj, zj, vxj, vyj,
     &         vzj, xjt, yjt, zjt, vxjt, vyjt, vzjt, istat, rstat, iud,
     &         trim(diro)//'/'//'discard.out', fopenstat, nleft)
         else
            nleft = ntp
         endif

c if it is time, output orb. elements,
         if (t .ge. tout) then

            if (btest(iflgchk,0).or.btest(iflgchk,1))  then ! bit 0 or 1 is set
               call io_write_frame_odea(t, nbod, ntp, oloc, matp, umatp,
     &            mat, umat, mass, eta, mu, xj, yj, zj, vxj, vyj, vzj,
     &            etatp, xjt, yjt, zjt, vxjt, vyjt, vzjt, istat,
     &            trim(diro)//'/'//outfile, iuf, fopenstat, iflgchk,
     &            matr)
            endif

            if (btest(iflgchk,2))  then ! bit 2 is set
               call anal_energy_write_odea(t, nbod, umat, mass, xj, yj,
     &            zj, vxj, vyj, vzj, iue, fopenstat, diro, eoff, energy,
     &            dt)
c               print*,'t',t,'DE/E=',(energy-energy0)/energy
            endif

            call io_oloc_write(t, nbod, ntp, oloc, oloct, iuo, iuotp,
     &           diro, fopenstat)

            tout = tout + dtout
         endif

c...     If it is time, do a dump
         if(t.ge.tdump) then

            tfrac = (t-t0)/(tstop-t0)
            write(*,998) t, tfrac, nleft
            call io_dump_pl_hjs(trim(diro)//'/'//'dump_pl.dat',
     &           nbod, oloc, mass, umat,
     &           xj, yj, zj, vxj, vyj, vzj, lclose, iflgchk, rplsq)
            call io_dump_tp_hjs(trim(diro)//'/'//'dump_tp.dat', nbod,
     &         ntp, matp, xjt, yjt, zjt, vxjt, vyjt, vzjt, istat, rstat)
            call io_dump_param(trim(diro)//'/'//'dump_param.dat',
     &           t, tstop, dt, dtout, dtdump, iflgchk,
     &           rmin, rmax, rmaxu, qmin, lclose, outfile)
            tdump = tdump + dtdump

         endif

      enddo
c------------ End of the big loop from time 't0' to time 'tstop'

c...  Do a final dump for possible resumption later

      call io_dump_pl_hjs(trim(diro)//'/'//'dump_pl.dat',
     &     nbod, oloc, mass, umat,
     &     xj, yj, zj, vxj, vyj, vzj, lclose, iflgchk, rplsq)
      call io_dump_tp_hjs(trim(diro)//'/'//'dump_tp.dat',
     &     nbod, ntp, matp,
     &     xjt, yjt, zjt, vxjt, vyjt, vzjt, istat, rstat)
      call io_dump_param(trim(diro)//'/'//'dump_param.dat',
     &     t, tstop, dt, dtout,
     &     dtdump, iflgchk, rmin, rmax, rmaxu, qmin, lclose, outfile)

      print*, 'DE/E=', (energy-energy0)/energy

      call util_exit(0)
      end ! odea
c-----------------------------------------------------------------------
