c***********************************************************************
c                        COORD_VB2VG_TP.F
c***********************************************************************
c Converts Barycentric coords to Generalized Jacobi
c (only one dimension)
c
c Inputs:
c    nbod ==> Number of bodies
c    matp  ==> Transformation matrix for tp
c    vxb, vyb, vzb ==> Barycentric particle velocities
c    vzbt, vybt, vzbt ==> Barycentric velocities of the tp
c
c Outputs:
c    vzjt, vyjt, vzjt ==> Generalized Jacobi velocity of the tp

      subroutine coord_vb2vg_tp(nbod, matp, vxb, vyb, vzb,
     &   vxbt, vybt, vzbt, vxjt, vyjt, vzjt)


      include '../swift.inc'

c...  Inputs:
      integer nbod
      real*8 matp(nbod)
      real*8 vxb(nbod), vyb(nbod), vzb(nbod)
      real*8 vxbt, vybt, vzbt

c...  Outputs:
      real*8 vxjt, vyjt, vzjt

c...  Internals:
      integer j
      real*8 sumvx, sumvy, sumvz

c----
c...  Executable code

c...  First calc. the array eta(*) then convert to jacobi coords

      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvz = 0.0d0

      do j = 1,nbod
         if (matp(j).ne.0.0d0) then
            sumvx = sumvx + matp(j)*vxb(j)
            sumvy = sumvy + matp(j)*vyb(j)
            sumvz = sumvz + matp(j)*vzb(j)
         end if
      end do

      vxjt = vxbt + sumvx
      vyjt = vybt + sumvy
      vzjt = vzbt + sumvz

      return
      end     ! coord_vb2vg_tp
c-----------------------------------------------------------------------
