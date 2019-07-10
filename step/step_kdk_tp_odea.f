c***********************************************************************
c                            STEP_KDK_TP_ODEA.F
c***********************************************************************
c This subroutine takes a step in Generalized Jacobi coords (HJS)
c Does a KICK then a DRIFT then a KICK.
c ONLY DOES TEST PARTICLES
c
c Input:
c    i1st                   ==>  = 0 if first step; = 1 not
c    nbod                   ==>  Number of massive bodies
c    ntp                    ==>  Number of test bodies
c    mass                   ==>  Masses of bodies
c    matp, umatp            ==>  Conversion vectors for tp's
c    oloct                  ==>  Link between tp's and orbits
c    eta, mu, etatp         ==>  Masses of centers & sats
c    xjbeg, yjbeg, zjbeg    ==>  Bodies Jac. position at beginning
c    xbbeg, ybbeg, zbbeg    ==>  Bodies bary position at beginning
c    vxjbeg, vyjbeg, vzjbeg ==>  Bodies Jac. veloc. at beginning
c    ir3jbeg, ir3jend       ==>  1/rj^3, beginning and end
c    axbbeg, aybbeg, azbbeg ==>  Bodies bary accs. at beginning
c    vxjh, vyjh, vzjh       ==>  Bodies velocities at middle point
c    xjend, yjend, zjend    ==>  Bodies Jac. position at end
c    xbend, ybend, zbend    ==>  Bodies bary position at end
c    axbend, aybend, azbend ==>  Bodies bary accs. at end
c    xjt, yjt, zjt          ==>  Initial tp position in Jac. coord
c    vxjt, vyjt, vzjt       ==>  Initial tp velocity in Jac. coord
c    istat                  ==>  Status of the test paricles
c             istat(i, 1) = 0     ==> active:  = 1 not
c             istat(i, 2) = -1    ==> Danby did not work
c    dt                     ==>  Time step
c
c Output:
c    xjt, yjt, zjt          ==>  Final position in Jac. coord
c    vxjt, vyjt, vzjt       ==>  Final position in Jac. coord

      subroutine step_kdk_tp_odea(i1st, nbod, ntp, matp, umatp, oloct,
     &   mass, eta, mu, etatp, xjbeg, yjbeg, zjbeg, vxjbeg, vyjbeg,
     &   vzjbeg, ir3jbeg, xbbeg, ybbeg, zbbeg, axbbeg, aybbeg, azbbeg,
     &   vxjh, vyjh, vzjh, xjend, yjend, zjend, ir3jend,
     &   xbend, ybend, zbend, axbend, aybend, azbend, xjt, yjt, zjt,
     &   vxjt, vyjt, vzjt, istat, dt, oloc, olocold, mat, umat, change,
     &   time)

      include '../swift.inc'

c...  Inputs Only:
      integer nbod, ntp, i1st, oloc(NPLMAX,NPLMAX), olocold(nbod, nbod)
      real*8 umat(NPLMAX, NTPMAX), mat(NPLMAX, NTPMAX)
      real*8 mass(nbod), eta(nbod), mu(nbod), dt, time
      real*8 xjbeg(nbod), yjbeg(nbod), zjbeg(nbod)
      real*8 vxjbeg(nbod), vyjbeg(nbod), vzjbeg(nbod)
      real*8 xbbeg(nbod), ybbeg(nbod), zbbeg(nbod)
      real*8 axbbeg(nbod), aybbeg(nbod), azbbeg(nbod)
      real*8 ir3jbeg(nbod), ir3jend(nbod), etatp(ntp)
      real*8 xjend(nbod), yjend(nbod), zjend(nbod)
      real*8 xbend(nbod), ybend(nbod), zbend(nbod)
      real*8 axbend(nbod), aybend(nbod), azbend(nbod)
      real*8 vxjh(nbod), vyjh(nbod), vzjh(nbod)
      logical change

c...  Inputs and Outputs:
      integer istat(NTPMAX, NSTAT), oloct(NPLMAX, NTPMAX)
      real*8 umatp(NPLMAX, NTPMAX), matp(NPLMAX, NTPMAX)
      real*8 xjt(ntp), yjt(ntp), zjt(ntp)
      real*8 vxjt(ntp), vyjt(ntp), vzjt(ntp)

c...  Internals:
      integer j, k, istati(1, NSTAT), i1sttp
      real*8 dth, axbttp, aybttp, azbttp, axjttp, ayjttp, azjttp
      real*8 axjt(NTPMAX), ayjt(NTPMAX), azjt(NTPMAX)
      real*8 xbttp, ybttp, zbttp, vxbttp, vybttp, vzbttp
      logical checkchange(NTPMAX)

      save axjt, ayjt, azjt, checkchange

c----
c...  Executable code

      dth = 0.5d0*dt

      i1sttp = 0

c...  Loop over all tp's
      do j = 1, ntp

         if (istat(j,1).eq.0) then

            if (i1st.eq.0) then

c...           Convert to barycentric coordinates
               call coord_g2b_tp(nbod, umatp(:,j), xjbeg, yjbeg, zjbeg,
     &            vxjbeg, vyjbeg, vzjbeg, xjt(j), yjt(j), zjt(j),
     &            vxjt(j), vyjt(j), vzjt(j), xbttp, ybttp, zbttp,
     &            vxbttp, vybttp, vzbttp)

c...           Get the acceleration in bary frame
               call getacch_tp_hjs(nbod, mass, eta, mu, xjbeg, yjbeg,
     &              zjbeg, xbbeg, ybbeg, zbbeg, ir3jbeg, oloct(:, j),
     &              etatp(j), xbttp, ybttp, zbttp, xjt(j), yjt(j),
     &              zjt(j), axbttp, aybttp, azbttp)

c...           Convert bary accels to Jacobi accels (use vel. procedure)
               call coord_vb2vg_tp(nbod, matp(:, j), axbbeg, aybbeg,
     &           azbbeg, axbttp, aybttp, azbttp, axjttp, ayjttp, azjttp)

               checkchange(j) = .False. ! Initial change

            else
               axjttp = axjt(j)
               ayjttp = ayjt(j)
               azjttp = azjt(j)
            end if

c...        Take into account a potential previous hierarchy change
c...        of massive bodies
            if (change) then

               call ce_change_hierarch_tp(i1sttp, nbod, mass, eta, mu,
     &            oloc, olocold, mat, umat, oloct(:, j), matp(:, j),
     &            umatp(:, j), etatp(j), xjbeg, yjbeg, zjbeg, vxjbeg,
     &            vyjbeg, vzjbeg, axbbeg, aybbeg, azbbeg, ir3jbeg,
     &            xjt(j), yjt(j), zjt(j), vxjt(j), vyjt(j), vzjt(j),
     &            axjttp, ayjttp, azjttp)

            end if

            if (checkchange(j)) then

               call ce_check_change_tp(j, nbod, mass, eta, mu, oloc,
     &              umat, oloct(:,j), matp(:,j), umatp(:,j), etatp(j),
     &              xjbeg, yjbeg, zjbeg, vxjbeg, vyjbeg, vzjbeg,
     &              axbbeg, aybbeg, azbbeg, ir3jbeg,
     &              xjt(j), yjt(j), zjt(j), vxjt(j), vyjt(j), vzjt(j),
     &              axjttp, ayjttp, azjttp, time)

            endif

c...        Apply a Jacobi kick for a half dt
            call kickvh_tp(1, vxjt(j), vyjt(j), vzjt(j),
     &         axjttp, ayjttp, azjttp, istat(j,1), dth)

c...        Take a drift forward full step
            do k=1,NSTAT
               istati(1, k) = istat(j, k)
            end do
            call drift_tp(1, etatp(j), xjt(j), yjt(j), zjt(j), vxjt(j),
     &         vyjt(j), vzjt(j), dt, istati)

c...        After drift, compute bary pos. and vels.
            call coord_g2b_tp(nbod, umatp(:,j), xjend, yjend, zjend,
     &         vxjh, vyjh, vzjh, xjt(j), yjt(j), zjt(j), vxjt(j),
     &         vyjt(j), vzjt(j), xbttp, ybttp, zbttp, vxbttp, vybttp,
     &         vzbttp)

c...        Get the acceleration in bary frame.
            call getacch_tp_hjs(nbod, mass, eta, mu, xjend, yjend,
     &         zjend, xbend, ybend, zbend, ir3jend, oloct(:,j),
     &         etatp(j), xbttp, ybttp, zbttp, xjt(j), yjt(j), zjt(j),
     &         axbttp, aybttp, azbttp)

c...        Convert bary accels to Jacobi accels (use vel. procedure)
            call coord_vb2vg_tp(nbod, matp(:,j), axbend, aybend, azbend,
     &         axbttp, aybttp, azbttp, axjttp, ayjttp, azjttp)
c...        Apply a Jacobi kick for a half dt
            call kickvh_tp(1, vxjt(j), vyjt(j), vzjt(j),
     &           axjttp, ayjttp, azjttp, istat(j,1), dth)

            call ce_crit_tp(xjt(j), yjt(j), zjt(j), etatp(j),
     &           axjttp, ayjttp, azjttp, checkchange(j))

            axjt(j) = axjttp
            ayjt(j) = ayjttp
            azjt(j) = azjttp

         end if

         i1sttp = 1

      end do

      i1st = 1

      return
      end   ! step_kdk_tp_odea
c-----------------------------------------------------------------------
