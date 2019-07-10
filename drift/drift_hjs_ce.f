c*************************************************************************
c                        DRIFT_HJS_CE.F
c*************************************************************************
c This subroutine loops through the particles and calls the danby routine
c    (Generalized Jacobi coordinates case)
c
c Input:
c    nbod        ==>  number of massive bodies (int scalar)
c    mass        ==>  mass of bodies (real array)
c    eta         ==> Masses of centers for orbits (real array)
c    mu          ==> Masses of satellites for orbits (real arr.)
c    xj,yj,zj    ==>  initial position in jacobi coord
c    vxj,vyj,vzj ==>  initial position in jacobi coord
c    dt          ==>  time step (real scalar)
c
c Output:
c    xj,yj,zj    ==>  final position in jacobi coord
c    vxj,vyj,vzj ==>  final position in jacobi coord

      subroutine drift_hjs_ce(nbod, mass, eta, mu, xj, yj, zj, vxj, vyj,
     &   vzj, dt, ce)

      include '../swift.inc'

c...  Inputs Only:
      integer nbod
      real*8 mass(nbod), eta(nbod), mu(nbod), dt
      logical ce(nbod)

c...  Inputs and Outputs:
      real*8 xj(nbod), yj(nbod), zj(nbod)
      real*8 vxj(nbod), vyj(nbod), vzj(nbod)

c...  Internals:
      real*8 gm
      integer j, iflg

c----
c...  Executable code

c...  Take a drift forward dth

      do j = 2, nbod
         if (.not.ce(j)) then
            gm = eta(j) + mu(j)
            call drift_one(gm, xj(j), yj(j), zj(j),
     &           vxj(j), vyj(j), vzj(j), dt, iflg)
            if(iflg.ne.0) then
               write(*,*) ' Orbit ',j,' is lost !'
               write(*,*) gm, dt
               write(*,*) xj(j), yj(j), zj(j)
               write(*,*) vxj(j), vyj(j), vzj(j)
               write(*,*) ' STOPPING '
               call util_exit(1)
            endif
         endif
      enddo

      return
      end ! drift_hjs_ce
c--------------------------------------------------------------------------
