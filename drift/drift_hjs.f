c***********************************************************************
c                        DRIFT_HJS.F
c***********************************************************************
c This subroutine compute the independent Keplerian motions
c  with the Danby routine (Generalized Jacobi coordinates case)
c
c Inputs:
c     nbod          ==> Number of massive bodies
c     mass          ==> Mass of bodies
c     eta           ==> Masses of centers for orbits
c     mu            ==> Masses of satellites for orbits
c     xj, yj ,zj    ==> Initial position in Jacobi coord
c     vxj, vyj, vzj ==> Initial velocity in Jacobi coord
c     dt            ==> Time step
c
c Outputs:
c    xj, yj, zj     ==>  Final position in Jacobi coord
c    vxj, vyj, vzj  ==>  Final velocity in Jacobi coord

      subroutine drift_hjs(nbod, mass, eta, mu, xj, yj, zj,
     &    vxj, vyj, vzj, dt)

      include '../swift.inc'

c...  Inputs Only:
      integer nbod
      real*8 mass(nbod), eta(nbod), mu(nbod), dt

c...  Inputs and Outputs:
      real*8 xj(nbod), yj(nbod), zj(nbod)
      real*8 vxj(nbod), vyj(nbod), vzj(nbod)

c...  Internals:
      real*8 gm
      integer j, iflg

c----
c...  Executable code

c...  Take a drift forward dt
      do j = 2,nbod
         gm = eta(j) + mu(j)
         call drift_one(gm, xj(j), yj(j), zj(j), vxj(j), vyj(j), vzj(j),
     &      dt, iflg)
         if(iflg.ne.0) then
            write(*,*) ' Orbit ',j,' is lost!'
            write(*,*) gm, dt
            write(*,*) xj(j), yj(j), zj(j)
            write(*,*) vxj(j), vyj(j), vzj(j)
            write(*,*) ' STOPPING '
            call util_exit(1)
         endif
      enddo

      return
      end ! drift_hjs
c-----------------------------------------------------------------------
