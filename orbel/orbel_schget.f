c***********************************************************************
c                      ORBEL_SCHGET.F
c***********************************************************************
c Purpose:  Given an angle, efficiently compute sinh and cosh.
c
c Input:
c   angle ==> angle in radians
c
c Outputs:
c   shx   ==>  sinh(angle)
c   chx   ==>  cosh(angle)
c
c Remarks: Based on the routine SCGET for sine's and cosine's.
c       We use the sqrt rather than cosh (it's faster)
c       Be sure the angle is in radians and it can't be larger than 300
c       or overflows will occur!

      subroutine orbel_schget(angle, shx, chx)

      include '../swift.inc'

c...  Inputs Only:
      real*8 angle

c...  Output:
      real*8 shx,chx

c----
c...  Executable code

      shx = sinh(angle)
      chx= sqrt(1.d0 + shx*shx)

      return
      end   ! orbel_schget
c-----------------------------------------------------------------------
