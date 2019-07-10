c***********************************************************************
c                            DRIFT_KEPMD.F
c***********************************************************************
c Subroutine for solving kepler's equation in difference form for an
c  ellipse, given SMALL dm and SMALL eccentricity.
c See DRIFT_DAN.F for the criteria.
c Built for speed : does not check how well the original
c  equation is solved! (can do that in the calling routine by
c  checking how close (x - ec*s +es*(1.-c) - dm) is to zero).
c
c Inputs:
c    dm      ==> Increment in mean anomaly M
c    es, ec  ==> Ecc. times sin and cos of E_0 (eccentric anomaly)
c
c Outputs:
c    x       ==> Solution to Kepler's difference equation
c    s, c    ==> Sine and cosine of x

      subroutine drift_kepmd(dm, es, ec, x, s, c)

      implicit none

c...  Inputs
      real*8 dm, es, ec

c...  Outputs
      real*8 x, s, c

c...  Internals
      real*8 A0, A1, A2, A3, A4
      parameter(A0 = 39916800.d0, A1 = 6652800.d0, A2 = 332640.d0)
      parameter(A3 = 7920.d0, A4 = 110.d0)
      real*8 dx
      real*8 fac1, fac2, q, y
      real*8 f, fp, fpp, fppp

c----
c...  Executable code

c...  Calc initial guess for root
      fac1 = 1.d0/(1.d0 - ec)
      q = fac1*dm
      fac2 = es*es*fac1 - ec/3.d0
      x = q*(1.d0 -0.5d0*fac1*q*(es -q*fac2))

c...  Excellent approx. to sin and cos of x for small x.
      y = x*x
      s = x*(A0-y*(A1-y*(A2-y*(A3-y*(A4-y)))))/A0
      c = sqrt(1.d0 - s*s)

c...  Compute better value for the root using quartic Newton method
      f = x - ec*s + es*(1.-c) - dm
      fp = 1. - ec*c + es*s
      fpp = ec*s + es*c
      fppp = ec*c - es*s
      dx = -f/fp
      dx = -f/(fp + 0.5*dx*fpp)
      dx = -f/(fp + 0.5*dx*fpp + 0.16666666666666666*dx*dx*fppp)
      x = x + dx

c...  Excellent approx. to sin and cos of x for small x.
      y = x*x
      s = x*(A0-y*(A1-y*(A2-y*(A3-y*(A4-y)))))/A0
      c = sqrt(1.d0 - s*s)

      return
      end ! drift_kepmd
c-----------------------------------------------------------------------
