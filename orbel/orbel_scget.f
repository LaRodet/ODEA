c***********************************************************************
c                            ORBEL_SCGET.F
c***********************************************************************
c Given an angle, efficiently compute sin and cos.
c Remarks: The HP 700 series won't return correct answers for sin
c  and cos if the angle is bigger than 3e7. We first reduce it
c  to the range [0,2pi) and use the sqrt rather than cos (it's faster)
c  BE SURE THE ANGLE IS IN RADIANS - NOT DEGREES!
c
c Input:
c    angle ==> angle in radians
c
c Outputs:
c    sx    ==>  sin(angle)
c    cx    ==>  cos(angle)

      subroutine orbel_scget(angle, sx, cx)

      include '../swift.inc'

c...  Inputs Only:
      real*8 angle

c...  Outputs:
      real*8 sx, cx

c... Internals:
      integer nper
      real*8 x
      real*8 PI3BY2
      parameter(PI3BY2 = 1.5d0*PI)

c----
c...  Executable code

      nper = angle/TWOPI
      x = angle - nper*TWOPI
      if(x.lt.0.d0) then
         x = x + TWOPI
      endif
      sx = sin(x)
      cx= sqrt(1.d0 - sx*sx)
      if( (x .gt. PIBY2) .and. (x .lt.PI3BY2)) then
         cx = -cx
      endif

      return
      end   ! orbel_scget
c-----------------------------------------------------------------------
