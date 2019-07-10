c*************************************************************************
c                            CE_CRIT_TP.F
c*************************************************************************
c This subroutine computes criteria to decide if the hierarchy of a tp
c should be checked
c
c Inputs:
c    etap          ==> Mass of centers of tp's orbit
c    xjt, yjt, zjt ==> Position of tp in Gen. Jacobi coord
c    axj, ayj, azj ==> Acceleration of tp in Gen. Jacobi coord
c
c Output:
c    checkchange   ==> True if a hierarchy check is needed

      subroutine ce_crit_tp(xjt, yjt, zjt, etatp, axjt, ayjt, azjt,
     &     checkchange)

      include '../swift.inc'
      include '../odea.inc'

c...  Inputs Only:
      real*8 etatp
      real*8 xjt, yjt, zjt
      real*8 axjt, ayjt, azjt

c...  Outputs Only
      logical checkchange

c...  Internals:
      real*8 r2, a, c, akep

c----
c...  Executable code

      checkchange = .False.

      r2 = xjt*xjt + yjt*yjt + zjt*zjt

      a = sqrt(axjt*axjt+ayjt*ayjt+azjt*azjt)
      akep = etatp/r2

      c = a/akep

      if (c.gt.CRITCHANGE) checkchange = .True.

      return
      end   ! ce_crit_c2
c-----------------------------------------------------------------------
