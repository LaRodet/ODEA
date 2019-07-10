c***********************************************************************
c                        DRIFT_KEPU_STUMPFF.F
c***********************************************************************
c Calculation of Stumpff functions
c see Danby p.172  equations 6.9.15
c
c Input:
c    x ==>  argument
c
c Outputs:
c    c0, c1, c2, c3   ==>  c's from p171-172

      subroutine drift_kepu_stumpff(x, c0, c1, c2, c3)

      include '../swift.inc'

c...  Inputs:
      real*8 x

c...  Outputs:
      real*8 c0, c1, c2, c3

c...  Internals:
      integer n, i
      real*8 xm

c----
c...  Executable code

      n = 0
      xm = 0.1d0
      do while(abs(x).ge.xm)
         n = n + 1
         x = x/4.0d0
      enddo

      c2 = (1.d0-x*(1.d0-x*(1.d0-x*(1.d0-x*(1.d0-x*(1.d0-x/182.d0)
     &       /132.d0)/90.d0)/56.d0)/30.d0)/12.d0)/2.d0
      c3 = (1.d0-x*(1.d0-x*(1.d0-x*(1.d0-x*(1.d0-x*(1.d0-x/210.d0)
     &       /156.d0)/110.d0)/72.d0)/42.d0)/20.d0)/6.d0
      c1 = 1.d0 - x*c3
      c0 = 1.d0 - x*c2

      if(n.ne.0) then
         do i=n,1,-1
            c3 = (c2 + c0*c3)/4.d0
            c2 = c1*c1/2.d0
            c1 = c0*c1
            c0 = 2.d0*c0*c0 - 1.d0
            x = x * 4.d0
         enddo
      endif

      return
      end     !   drift_kepu_stumpff
c-----------------------------------------------------------------------
