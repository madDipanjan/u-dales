!> \file modboundary.f90
!! All boundary conditions are in this file, except for immersed boundaries.
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!
!! \par Revision list
!! \par Authors
!!
module modboundary

   implicit none
   save
   private
   public :: initboundary, boundary, grwdamp, ksp, tqaver, halos, bcp, bcpup, closurebc, &
             xm_periodic, xT_periodic, xq_periodic, xs_periodic, ym_periodic, yT_periodic, yq_periodic, ys_periodic
   integer :: ksp = -1 !<    lowest level of sponge layer
   real, allocatable :: tsc(:) !<   damping coefficients to be used in grwdamp.
   real :: rnu0 = 2.75e-3
contains
   !>
   !! Initializing Boundary; specifically the sponge layer
   !>
   subroutine initboundary
      use modglobal,    only : ib, kb, ke, kh, kmax, pi, zf, iplane
      use modinletdata, only : irecy
      implicit none

      real    :: zspb, zspt
      integer :: k
      allocate (tsc(kb:ke + kh))
      ! Sponge layer
      if (ksp == -1) then
         !      ksp  = min(3*kmax/4,kmax - 15)
         ksp = (kb - 1) + max(min(3*kmax/4, kmax - 15),1)
      end if

      zspb = zf(ksp)
      zspt = zf(ke)

      tsc(kb:ksp - 1) = 0.0
      do k = ksp, ke
         tsc(k) = rnu0*sin(0.5*pi*(zf(k) - zspb)/(zspt - zspb))**2
      end do
      do k = 1, kh
        tsc(ke + k) = tsc(ke + k - 1)
      end do
      irecy = ib + iplane

   end subroutine initboundary

   !>
   !! Fill halo cells, including ghost cells outside domain
   ! Needs to be called before divergence is calculated
   subroutine halos

      use modglobal, only : ib, ie, ih, jb, je, jh, kb, ke, kh, ihc, jhc, khc, nsv, &
                            BCxm, BCym, BCxT, BCyT, BCxq, BCyq, BCxs, BCys, &
                            BCxm_periodic, BCxT_periodic, BCxq_periodic, BCxs_periodic, &
                            BCym_periodic, BCyT_periodic, BCyq_periodic, BCys_periodic, &
                            ibrank, ierank, jbrank, jerank
      use modfields, only : u0, v0, w0, um, vm, wm, thl0, thlm, qt0, qtm, sv0, svm, thl0c
      use decomp_2d, only : exchange_halo_z
      implicit none
      integer i, k, n

      call exchange_halo_z(u0)
      call exchange_halo_z(v0)
      call exchange_halo_z(w0)
      call exchange_halo_z(um)
      call exchange_halo_z(vm)
      call exchange_halo_z(wm)
      call exchange_halo_z(thl0)
      call exchange_halo_z(thlm)
      call exchange_halo_z(thl0c, opt_zlevel=(/ihc,jhc,khc/))
      call exchange_halo_z(qt0)
      call exchange_halo_z(qtm)
      do n = 1, nsv
         call exchange_halo_z(sv0(:, :, :, n), opt_zlevel=(/ihc,jhc,khc/))
         call exchange_halo_z(svm(:, :, :, n), opt_zlevel=(/ihc,jhc,khc/))
      enddo

      if (ibrank .and. ierank) then ! not parallelized in x
        if (BCxm == BCxm_periodic) call xm_periodic
        if (BCxT == BCxT_periodic) call xT_periodic
        if (BCxq == BCxq_periodic) call xq_periodic
        if (BCxs == BCxs_periodic) call xs_periodic
      end if

      if (jbrank .and. jerank) then ! not parallelized in x
        if (BCym == BCym_periodic) call ym_periodic
        if (BCyT == BCyT_periodic) call yT_periodic
        if (BCyq == BCyq_periodic) call yq_periodic
        if (BCys == BCys_periodic) call ys_periodic
      end if

    end subroutine halos


   !>
   !! Set boundary conditions for the next timestep
   ! Will result in velocity field being not divergence-free
   subroutine boundary
      use modglobal,      only : ib, ie, ih, jb, je, jh, kb, ke, kh, ihc, jhc, khc, dzf, zh, nsv, &
                                 ltempeq, lmoist, luvolflowr, luoutflowr, &
                                 BCxm, BCym, BCxT, BCyT, BCxq, BCyq, BCxs, BCys, BCtopm, BCtopT, BCtopq, BCtops, &
                                 BCtopm_freeslip, BCtopm_noslip, BCtopm_pressure, &
                                 BCtopT_flux, BCtopT_value, BCtopq_flux, BCtopq_value, BCtops_flux, BCtops_value, &
                                 BCxm_periodic, BCxm_profile, BCxm_driver, &
                                 BCxT_periodic, BCxT_profile, BCxT_driver, &
                                 BCxq_periodic, BCxq_profile, BCxq_driver, &
                                 BCxs_periodic, BCxs_profile, BCxs_driver, BCxs_custom, &
                                 BCym_periodic, BCym_profile, BCyT_periodic, BCyT_profile, &
                                 BCyq_periodic, BCyq_profile, BCys_periodic, &
                                 ibrank, ierank, jbrank, jerank, e12min, idriver, &
                                 Uinf, Vinf, &
                                 rk3step, lchunkread
      use modfields,      only : u0, v0, w0, um, vm, wm, thl0, thlm, qt0, qtm, e120, e12m, sv0, svm, u0av, v0av, uouttot, vouttot, thl0c
      use modsubgriddata, only : ekh, ekm, loneeqn
      use modsurfdata,    only : thl_top, qt_top, sv_top, wttop, wqtop, wsvtop
      use modmpi,         only : myid, slabsum, avey_ibm
      use moddriver,      only : drivergen, driverchunkread
      use modinletdata,   only : ubulk, vbulk, iangle
      use decomp_2d,      only : exchange_halo_z

      implicit none
      real, dimension(kb:ke) :: uaverage, vaverage
      real, dimension(ib:ie,kb:ke) :: uavey
      integer i, k, n

     ! if not using massflowrate need to set outflow velocity
     if (luoutflowr) then
        ! do nothing - calculated in modforces
     elseif (.not. luvolflowr) then
        !ubulk = sum(u0av)/(ke-kb+1)
        do k = kb, ke
           uaverage(k) = u0av(k)*dzf(k)
        end do

        do k = kb, ke
           vaverage(k) = v0av(k)*dzf(k)
        end do
        ! need a method to know if we have all blocks at lowest cell kb
        ! assuming this for now (hence kb+1)
        uouttot = sum(uaverage(kb:ke))/(zh(ke + 1) - zh(kb))
        vouttot = sum(vaverage(kb:ke))/(zh(ke + 1) - zh(kb))
     else
        uouttot = ubulk
        vouttot = vbulk
     end if

     ! Bottom BC - many ways of enforcing this but this is simplest
     ! Other variables handled by bottom
     wm(:, :, kb) = 0.
     w0(:, :, kb) = 0.

     !! Top
     ! Momentum
     select case(BCtopm)
     case(BCtopm_freeslip)
        !free-slip = zero-flux
        call fluxtop(um, ekm, 0.0)
        call fluxtop(u0, ekm, 0.0)
        call fluxtop(vm, ekm, 0.0)
        call fluxtop(v0, ekm, 0.0)
        do k = 1, kh
          w0(:, :, ke + k) = 0.0
          wm(:, :, ke + k) = 0.0
          if (loneeqn) then
            e120(:, :, ke + k) = e12min
            e12m(:, :, ke + k) = e12min
          end if
        end do
     case(BCtopm_noslip)
        !no-slip = fixed velocity at wall
        call valuetop(um, Uinf)
        call valuetop(u0, Uinf)
        call valuetop(vm, Vinf)
        call valuetop(v0, Vinf)
        do k = 1, kh
          w0(:, :, ke + k) = 0.0
          wm(:, :, ke + k) = 0.0
        end do
      case(BCtopm_pressure)
         call fluxtop(um, ekm, 0.0)
         call fluxtop(u0, ekm, 0.0)
         call fluxtop(vm, ekm, 0.0)
         call fluxtop(v0, ekm, 0.0)
         if (loneeqn) then
          do k = 1, kh
           e120(:, :, ke + k) = e12min
           e12m(:, :, ke + k) = e12min
          end do
         end if
         ! w considered in modpois
      case default
        write(0, *) "ERROR: top boundary type for velocity undefined"
        stop 1
     end select

     ! Temperature
     select case(BCtopT)
     case(BCtopT_flux)
        call fluxtop(thlm, ekh, wttop)
        call fluxtop(thl0, ekh, wttop)
        do n=1,khc
           thl0c(:,:,ke+n) = thl0c(:,:,ke+n-1)
        end do
     case(BCtopT_value)
        call valuetop(thlm, thl_top)
        call valuetop(thl0, thl_top)
     case default
        write(0, *) "ERROR: top boundary type for temperature undefined"
        stop 1

     end select

     ! Moisture
     select case(BCtopq)
     case(BCtopq_flux)
        call fluxtop(qtm, ekh, wqtop)
        call fluxtop(qt0, ekh, wqtop)
     case(BCtopq_value)
        call valuetop(qtm, qt_top)
        call valuetop(qt0, qt_top)
     case default
        write(0, *) "ERROR: top boundary type for moisture undefined"
        stop 1
     end select

     ! Scalars
     select case(BCtops)
     case(BCtops_flux)
        call fluxtopscal(wsvtop)
        call fluxtopscal(wsvtop)
     case(BCtops_value)
        call valuetopscal(sv_top)
        call valuetopscal(sv_top)
     case default
        write(0, *) "ERROR: top boundary type for scalars undefined"
        stop 1
     end select

     if (idriver == 1) call drivergen ! Should be moved elsewhere, as not related to boundary conditions.

     ! x inlet
     if (ibrank) then ! set inlet
       ! Momentum
       select case(BCxm)
       case(BCxm_periodic)
         ! Handled in halos
       case(BCxm_profile)
         !uouttot = cos(iangle)*ubulk
         call xmi_profile
       case(BCxm_driver)
         !uouttot = ubulk ! does this hold for all forcings of precursor simulations? tg3315
         if(rk3step==0 .or. rk3step==3) then
          if (lchunkread) call driverchunkread
          call drivergen ! think this should be done at the start of an rk3 loop?
         end if
         call xmi_driver
       case default
         write(0, *) "ERROR: lateral boundary type for veloctiy in x-direction undefined"
         stop 1
       end select

       ! Temperature
       if (ltempeq) then
         select case(BCxT)
         case(BCxT_periodic) ! periodic
           ! Handled in halos
         case(BCxT_profile) ! profile
           call xTi_profile
         case(BCxT_driver)
           call xTi_driver
         case default
           write(0, *) "ERROR: lateral boundary type for temperature in x-direction undefined"
           stop 1
         end select
       end if

       ! Moisture
       if (lmoist) then
         select case(BCxq)
         case(BCxq_periodic)
           ! Handled in halos
         case(BCxq_profile)
           call xqi_profile
         case(BCxq_driver)
           call xqi_driver
         case default
           write(0, *) "ERROR: lateral boundary type for humidity in x-direction undefined"
           stop 1
         end select
       end if

       ! Scalars
       if (nsv > 0) then
         select case(BCxs)
         case(BCxs_periodic)
           ! Handled in halos
         case(BCxs_profile)
           call xsi_profile
         case(BCxs_driver)
           call xsi_driver
        case(BCxs_custom)
           call xsi_custom
         case default
           write(0, *) "ERROR: lateral boundary type for scalars in x-direction undefined"
           stop 1
         end select
       end if

     end if !ibrank

     if (jbrank) then ! set y inlet
       ! Momentum
       select case(BCym)
       case(BCym_periodic)
         ! Handled in halos
       case(BCym_profile)
         call ymi_profile
       case default
         write(0, *) "ERROR: lateral boundary type for veloctiy in y-direction undefined"
         stop 1
       end select

       ! Temperature
       if (ltempeq) then
         select case(BCyT)
         case(BCyT_periodic)
           ! Handled in halos
         case(BCyT_profile)
           call yTi_profile
         case default
           write(0, *) "ERROR: lateral boundary type for temperature in y-direction undefined"
           stop 1
         end select
       end if

       ! Moisture
       if (lmoist) then
         select case(BCyq)
         case(BCyq_periodic)
           ! Handled in halos
         case(BCyq_profile)
           call yqi_profile
         case default
           write(0, *) "ERROR: lateral boundary type for humidity in y-direction undefined"
           stop 1
         end select
       end if

       if (nsv > 0) then !scalars
         select case(BCys)
         case(1)
           ! Handled in halos
         case(2)
           call ysi_profile
         case default
           write(0, *) "ERROR: lateral boundary type for scalars in y-direction undefined"
           stop 1
         end select
       end if

     end if !jbrank

     !> Outlet
     ! Currently only outflow boundary conditions are convective
     if (ierank) then
       if (BCxm .ne. BCxm_periodic) call xmo_convective
       if ((BCxT .ne. BCxT_periodic) .and. ltempeq) call xTo_convective
       if ((BCxq .ne. BCxq_periodic) .and. lmoist ) call xqo_convective
       if ((BCxs .ne. BCxs_periodic) .and. nsv > 0) call xso_convective
     end if

     if (jerank) then
       if (BCym .ne. BCym_periodic) call ymo_convective
       if ((BCyT .ne. BCyT_periodic) .and. ltempeq) call yTo_convective
       if ((BCyq .ne. BCyq_periodic) .and. lmoist ) call yqo_convective
       if ((BCys .ne. BCys_periodic) .and. nsv > 0) call yso_convective
     end if

   end subroutine boundary


   subroutine closurebc
     use modsubgriddata, only : ekm, ekh
     use modglobal,      only : ib, ie, jb, je, kb, ke, ih, jh, kh, numol, prandtlmoli, &
                                ibrank, ierank, jbrank, jerank, BCtopm, BCxm, BCym, &
                                BCtopm_freeslip, BCtopm_noslip, BCtopm_pressure, &
                                BCxm_periodic, BCym_periodic
     use decomp_2d,      only : exchange_halo_z
     integer i, j, k

     call exchange_halo_z(ekm)
     call exchange_halo_z(ekh)

     ! Top and bottom
     if ((BCtopm .eq. BCtopm_freeslip) .or. (BCtopm .eq. BCtopm_pressure)) then
       do j = jb - jh, je + jh
         do i = ib - ih, ie + ih
          do k = 1, kh
           ekm(i, j, ke + k) = ekm(i, j, ke + k - 1) ! zero-gradient top wall
           ekh(i, j, ke + k) = ekh(i, j, ke + k - 1) ! zero-gradient top wall
           ekm(i, j, kb - k) = 2.*numol - ekm(i, j, kb - k + 1) ! no-slip lower wall
           ekh(i, j, kb - k) = (2.*numol*prandtlmoli) - ekh(i, j, kb - k + 1) ! no-slip lower wall
          end do
         end do
       end do
     else if (BCtopm .eq. BCtopm_noslip) then
       do j = jb - jh, je + jh
         do i = ib - ih, ie + ih
          do k = 1, kh
           ekm(i, j, ke + k) = 2.*numol - ekm(i, j, ke + k - 1) ! no-slip top wall
           ekh(i, j, ke + k) = (2.*numol*prandtlmoli) - ekh(i, j, ke + k - 1) ! no-slip top wall
           ekm(i, j, kb - k) = 2.*numol - ekm(i, j, kb - k + 1) ! no-slip lower wall
           ekh(i, j, kb - k) = (2.*numol*prandtlmoli) - ekh(i, j, kb - k + 1) ! no-slip lower wall
          end do
         end do
       end do
     end if

     if (BCxm .ne. BCxm_periodic) then ! inflow/outflow
       if (ibrank) then
        do i = 1, ih
         ekm(ib - i, :, :) = ekm(ib - i + 1, :, :)
         ekh(ib - i, :, :) = ekh(ib - i + 1, :, :)
        end do
       end if
       if (ierank) then
        do i = 1, ih
         ekm(ie + i, :, :) = ekm(ie + i - 1, :, :)
         ekh(ie + i, :, :) = ekh(ie + i - 1, :, :)
        end do
       end if
     else ! periodic
       if (ibrank .and. ierank) then
        do i = 1, ih
         ekm(ib - i, :, :) = ekm(ie + i - 1, :, :)
         ekm(ie + i, :, :) = ekm(ib - i + 1, :, :)
         ekh(ib - i, :, :) = ekh(ie + i - 1, :, :)
         ekh(ie + i, :, :) = ekh(ib - i + 1, :, :)
        end do
       end if
     end if

     if (BCym .ne. BCym_periodic) then ! inflow/outflow
       if (jbrank) then
        do j = 1, jh
         ekm(:,jb-j,:) = ekm(:,jb-j+1,:)
         ekh(:,jb-j,:) = ekh(:,jb-j+1,:)
        end do
       end if
       if (jerank) then
        do j = 1, jh
         ekm(:,je+j,:) = ekm(:,je+j-1,:)
         ekh(:,je+j,:) = ekh(:,je+j-1,:)
        end do
       end if
     else ! periodic
       if (jbrank .and. jerank) then
        do j = 1, jh
         ekm(:, jb - j, :) = ekm(:, je + j - 1, :)
         ekm(:, je + j, :) = ekm(:, jb - j + 1, :)
         ekh(:, jb - j, :) = ekh(:, je + j - 1, :)
         ekh(:, je + j, :) = ekh(:, jb - j + 1, :)
        end do
       end if
     end if

   end subroutine closurebc

   !>set lateral periodic boundary conditions for momentum in x/i direction
   subroutine xm_periodic
      use modglobal, only : ib, ie, ih
      use modfields, only : u0, um, v0, vm, w0, wm, e120, e12m
      use modsubgriddata, only : loneeqn, lsmagorinsky
      use modmpi, only : excis

      integer n, m

      do m = 1, ih
         u0(ib - m, :, :) = u0(ie + 1 - m, :, :)
         u0(ie + m, :, :) = u0(ib - 1 + m, :, :)
         v0(ib - m, :, :) = v0(ie + 1 - m, :, :)
         v0(ie + m, :, :) = v0(ib - 1 + m, :, :)
         w0(ib - m, :, :) = w0(ie + 1 - m, :, :)
         w0(ie + m, :, :) = w0(ib - 1 + m, :, :)
         um(ib - m, :, :) = um(ie + 1 - m, :, :)
         um(ie + m, :, :) = um(ib - 1 + m, :, :)
         vm(ib - m, :, :) = vm(ie + 1 - m, :, :)
         vm(ie + m, :, :) = vm(ib - 1 + m, :, :)
         wm(ib - m, :, :) = wm(ie + 1 - m, :, :)
         wm(ie + m, :, :) = wm(ib - 1 + m, :, :)

         if (loneeqn) then
          e120(ib - m, :, :) = e120(ie + 1 - m, :, :)
          e120(ie + m, :, :) = e120(ib - 1 + m, :, :)
          e12m(ib - m, :, :) = e12m(ie + 1 - m, :, :)
          e12m(ie + m, :, :) = e12m(ib - 1 + m, :, :)
         end if
      end do

      return
   end subroutine xm_periodic


   !> Sets x/i periodic boundary conditions for the temperature
   subroutine xT_periodic
      use modglobal, only : ib, ie, ih, ihc
      use modfields, only : thl0, thlm, thl0c
      integer m

      do m = 1, ih
         thl0(ib - m, :, :) = thl0(ie + 1 - m, :, :)
         thl0(ie + m, :, :) = thl0(ib - 1 + m, :, :)
         thlm(ib - m, :, :) = thlm(ie + 1 - m, :, :)
         thlm(ie + m, :, :) = thlm(ib - 1 + m, :, :)
      end do

      do m = 1, ihc
         thl0c(ib - m, :, :) = thl0c(ie + 1 - m, :, :)
         thl0c(ie + m, :, :) = thl0c(ib - 1 + m, :, :)
      end do

      return
   end subroutine xT_periodic

   !> Sets x/i periodic boundary conditions for the humidity
   subroutine xq_periodic
      use modglobal, only : ib, ie, ih
      use modfields, only : qt0, qtm
      integer m

      do m = 1, ih
         qt0(ib - m, :, :) = qt0(ie + 1 - m, :, :)
         qt0(ie + m, :, :) = qt0(ib - 1 + m, :, :)
         qtm(ib - m, :, :) = qtm(ie + 1 - m, :, :)
         qtm(ie + m, :, :) = qtm(ib - 1 + m, :, :)
      end do

      return
   end subroutine xq_periodic

   !> Sets x/iperiodic boundary conditions for the scalars
   subroutine xs_periodic
      use modglobal, only : ib, ie, ihc
      use modfields, only : sv0, svm
      integer m, n

      do m = 1, ihc
         sv0(ib - m, :, :, :) = sv0(ie + 1 - m, :, :, :)
         sv0(ie + m, :, :, :) = sv0(ib - 1 + m, :, :, :)
         svm(ib - m, :, :, :) = svm(ie + 1 - m, :, :, :)
         svm(ie + m, :, :, :) = svm(ib - 1 + m, :, :, :)
      end do

      return
   end subroutine xs_periodic

   !>set lateral periodic boundary conditions for momentum in y/j direction
   subroutine ym_periodic
      use modglobal, only:ib, ie, jb, je, ih, jh, kb, ke, kh, jmax
      use modfields, only:u0, um, v0, vm, w0, wm, e120, e12m, shear
      use modsubgriddata, only:loneeqn, lsmagorinsky
      use modmpi, only:excjs

      integer n, m

      do m = 1, jh
         u0(:, jb - m, :) = u0(:, je + 1 - m, :)
         u0(:, je + m, :) = u0(:, jb - 1 + m, :)
         v0(:, jb - m, :) = v0(:, je + 1 - m, :)
         v0(:, je + m, :) = v0(:, jb - 1 + m, :)
         w0(:, jb - m, :) = w0(:, je + 1 - m, :)
         w0(:, je + m, :) = w0(:, jb - 1 + m, :)
         um(:, jb - m, :) = um(:, je + 1 - m, :)
         um(:, je + m, :) = um(:, jb - 1 + m, :)
         vm(:, jb - m, :) = vm(:, je + 1 - m, :)
         vm(:, je + m, :) = vm(:, jb - 1 + m, :)
         wm(:, jb - m, :) = wm(:, je + 1 - m, :)
         wm(:, je + m, :) = wm(:, jb - 1 + m, :)

         if (loneeqn) then
          e120(:, jb - m, :) = e120(:, je + 1 - m, :)
          e120(:, je + m, :) = e120(:, jb - 1 + m, :)
          e12m(:, jb - m, :) = e12m(:, je + 1 - m, :)
          e12m(:, je + m, :) = e12m(:, jb - 1 + m, :)
        end if
      end do

      return
   end subroutine ym_periodic


   !> Sets y/j periodic boundary conditions for the temperature
   subroutine yT_periodic
      use modglobal, only : jb, je, jh, jhc
      use modfields, only : thl0, thlm, thl0c
      use modmpi, only:excjs, myid, nprocs
      integer m

      do m = 1, jh
         thl0(:, jb - m, :) = thl0(:, je + 1 - m, :)
         thl0(:, je + m, :) = thl0(:, jb - 1 + m, :)
         thlm(:, jb - m, :) = thlm(:, je + 1 - m, :)
         thlm(:, je + m, :) = thlm(:, jb - 1 + m, :)
      end do

      do m = 1, jhc
         thl0c(:, jb - m, :) = thl0c(:, je + 1 - m, :)
         thl0c(:, je + m, :) = thl0c(:, jb - 1 + m, :)
      end do

      return
   end subroutine yT_periodic

   !> Sets y/j periodic boundary conditions for the humidity
   subroutine yq_periodic
      use modglobal, only : jb, je, jh
      use modfields, only : qt0, qtm

      integer m

      do m = 1, jh
         qt0(:, jb - m, :) = qt0(:, je + 1 - m, :)
         qt0(:, je + m, :) = qt0(:, jb - 1 + m, :)
         qtm(:, jb - m, :) = qtm(:, je + 1 - m, :)
         qtm(:, je + m, :) = qtm(:, jb - 1 + m, :)
      end do

      return
   end subroutine yq_periodic

   !> Sets y/j periodic boundary conditions for the scalars
   subroutine ys_periodic
      use modglobal, only : jb, je, jhc, nsv
      use modfields, only : sv0, svm
      integer n, m

      do n = 1, nsv
        do m = 1, jhc
          sv0(:, jb - m, :, :) = sv0(:, je + 1 - m, :, :)
          sv0(:, je + m, :, :) = sv0(:, jb - 1 + m, :, :)
          svm(:, jb - m, :, :) = svm(:, je + 1 - m, :, :)
          svm(:, je + m, :, :) = svm(:, jb - 1 + m, :, :)
        end do
      end do

      return
   end subroutine ys_periodic


     subroutine xmi_profile
       use modglobal,      only : ib, ie, jb, je, kb, ke, ih, jh, kh
       use modfields,      only : u0, um, v0, vm, w0, wm, e120, e12m, uprof, vprof, e12prof
       use modsubgriddata, only : loneeqn

       integer i, j, k

       do j = jb - jh, je + jh
         do k = kb, ke + kh
          u0(ib, j, k) = uprof(k)
          um(ib, j, k) = uprof(k)
          do i = 1, ih
           u0(ib - i, j, k) = 2*u0(ib-i+1, j, k) - u0(ib-i+2, j, k) ! (u(ib+1)+u(ib-1))/2 = u(ib)
           um(ib - i, j, k) = 2*um(ib-i+1, j, k) - um(ib-i+2, j, k) ! (u(ib+1)+u(ib-1))/2 = u(ib)
           if (i==1) then
            v0(ib - i, j, k) = 2*vprof(k) - v0(ib-i+1, j, k) ! (v(ib)+v(ib-1))/2 = vprof
            vm(ib - i, j, k) = 2*vprof(k) - vm(ib-i+1, j, k) ! (v(ib)+v(ib-1))/2 = vprof
           else
            v0(ib - i, j, k) = 2*v0(ib-i+1, j, k) - v0(ib-i+2, j, k) ! (v(ib-2)+v(ib))/2 = v(ib-1)
            vm(ib - i, j, k) = 2*vm(ib-i+1, j, k) - vm(ib-i+2, j, k) ! (v(ib-2)+v(ib))/2 = v(ib-1)
           end if
           w0(ib - i, j, k) = -w0(ib-i+1, j, k)
           wm(ib - i, j, k) = -wm(ib-i+1, j, k)
          end do
         end do
       end do

       if (loneeqn) then
         do j = jb - jh, je + jh
           do k = kb, ke + kh
            do i = 1, ih
              if (i==1) then
                e120(ib - i, j, k) = 2*e12prof(k) - e120(ib - i + 1, j, k) ! (e12(ib)+e12(ib-1))/2=e12prof
                e12m(ib - i, j, k) = 2*e12prof(k) - e12m(ib - i + 1, j, k) ! (e12(ib)+e12(ib-1))/2=e12prof
              else
                e120(ib - i, j, k) = 2*e120(ib-i+1, j, k) - e120(ib-i+2, j, k) ! (e12(ib)+e12(ib-2))/2=e12(ib-1)
                e12m(ib - i, j, k) = 2*e12m(ib-i+1, j, k) - e12m(ib-i+2, j, k) ! (e12(ib)+e12(ib-2))/2=e12(ib-1)
              endif
            end do
           end do
         end do
       end if

     end subroutine xmi_profile


     subroutine xmi_driver
       use modglobal,      only : ib, ie, jb, je, kb, ke, ih, jh, kh
       use modinletdata,   only : u0driver, umdriver, v0driver, vmdriver, w0driver, wmdriver
       use modfields,      only : u0, um, v0, vm, w0, wm, e120, e12m, e12prof
       use modsubgriddata, only : loneeqn

       integer i, j, k

        do i = 0,ih
          do j = jb - jh, je + jh
            do k = kb, ke !tg3315 removed +1 following above...
              u0(ib-i,j,k) = u0driver(j,k) !max(0.,u0driver(j,k))
              um(ib-i,j,k) = umdriver(j,k) !max(0.,umdriver(j,k))
              ! u0(ib,j,k) = u0driver(j,k) !max(0.,u0driver(j,k))
              ! um(ib,j,k) = umdriver(j,k) !max(0.,umdriver(j,k))
              ! u0(ib-1,j,k) = u0driver(j,k) !max(0.,u0driver(j,k))
              ! um(ib-1,j,k) = umdriver(j,k) !max(0.,umdriver(j,k))
              !! u0(ib-1,j,k)= 2*u0(ib, j, k) - u0(ib + 1, j, k) ! (u(ib+1)+u(ib-1))/2 = u(ib)
              !! um(ib-1,j,k)= 2*um(ib, j, k) - um(ib + 1, j, k) ! (u(ib+1)+u(ib-1))/2 = u(ib)
            end do
          end do
        end do

        do i = 1,ih
          do j = jb - jh, je + jh
            do k = kb, ke !tg3315 removed +1 following above...
              v0(ib-i,j,k)   = v0driver(j,k) !max(0.,v0driver(j,k))
              vm(ib-i,j,k)   = vmdriver(j,k) !max(0.,vmdriver(j,k))
            end do

            do k=kb,ke+kh
              w0(ib-i,j,k) = w0driver(j,k) !max(0.,w0driver(j,k))
              wm(ib-i,j,k) = wmdriver(j,k) !max(0.,wmdriver(j,k))
            end do
          end do
        end do

       if (loneeqn) then
        do i = 1,ih
         do j = jb - jh, je + jh
           do k = kb, ke + kh
             ! to be changed in the future: e12 should be taken from recycle plane!
             !e120(ib-1,j,k) = e120driver(j,k)      ! extrapolate e12 from interior
             !e12m(ib-1,j,k) = e12mdriver(j,k)      ! extrapolate e12 from interior
             e120(ib - i, j, k) = e120(ib-i+1, j, k) ! (e12(ib)+e12(ib-1))/2=e12prof
             e12m(ib - i, j, k) = e12m(ib-i+1, j, k) ! (e12(ib)+e12(ib-1))/2=e12prof
           end do
         end do
        end do
       end if

     end subroutine xmi_driver


     subroutine xTi_profile
       use modglobal, only : ib, ie, jb, je, kb, ke, ih, jh, kh
       use modfields, only : thl0, thlm, thlprof
       integer i, j, k

       ! set ghost cell
       ! do j = jb - 1, je + 1
       !   do k = kb, ke + 1
       !     thl0(ib - 1, j, k) = 2*thlprof(k) - thl0(ib, j, k)
       !     thlm(ib - 1, j, k) = 2*thlprof(k) - thlm(ib, j, k)
       !   end do
       ! end do
       do j = jb - jh, je + jh
         do k = kb, ke + kh
          do i = 0, ih
            thl0(ib - i, j, k) = thlprof(k)
            thlm(ib - i, j, k) = thlprof(k)
          end do
         end do
       end do

      !  ! set first internal cell as well
      !  do j = jb - jh, je + jh
      !   do k = kb, ke
      !      thl0(ib, j, k) = thlprof(k)
      !      thlm(ib, j, k) = thlprof(k)
      !   end do
      !  end do

     end subroutine xTi_profile


     subroutine xTi_driver
       use modglobal,    only : ib, ie, jb, je, kb, ke, ih, jh, kh
       use modinletdata, only : thl0driver, thlmdriver
       use modfields,    only : thl0, thlm
       integer i, j, k

       do j = jb - jh, je + jh
         do k = kb, ke + kh
          do i = 1, ih
           thl0(ib - i, j, k) = thl0driver(j, k)
           thlm(ib - i, j, k) = thlmdriver(j, k)
          end do
         end do
       end do

     end subroutine xTi_driver


     subroutine xqi_profile
       use modglobal,    only : ib, ie, jb, je, kb, ke, ih, jh, kh
       use modfields,    only : qt0, qtm, qtprof
       integer i, j, k

       do j = jb - jh, je + jh
         do k = kb, ke + kh
          do i = 1, ih
            if (i==1) then
              qt0(ib - i, j, k) = 2*qtprof(k) - qt0(ib - i + 1, j, k)
              qtm(ib - i, j, k) = 2*qtprof(k) - qtm(ib - i + 1, j, k)
            else
              qt0(ib - i, j, k) = 2*qt0(ib - i + 1, j, k) - qt0(ib - i + 2, j, k)
              qtm(ib - i, j, k) = 2*qtm(ib - i + 1, j, k) - qtm(ib - i + 2, j, k)
            end if
          end do
         end do
       end do

   end subroutine xqi_profile


   subroutine xqi_driver
     use modglobal,    only : ib, ie, jb, je, kb, ke, ih, jh, kh
     use modinletdata, only : qt0driver, qtmdriver
     use modfields,    only : qt0, qtm

     integer i, j, k

     do j = jb - jh, je + jh
       do k = kb, ke + kh
        do i = 1, ih
         qt0(ib - i, j, k) = qt0driver(j, k)
         qtm(ib - i, j, k) = qtmdriver(j, k)
        end do
       end do
     end do

   end subroutine xqi_driver


   subroutine xsi_profile
     use modglobal,    only : ib, ie, jb, je, kb, ke, nsv, ihc, khc
     use modfields,    only : sv0, svm, svprof

     integer j, k, n, m

     do j = jb, je
       do k = kb, ke + khc
         do n = 1, nsv
           do m = 1, ihc
             sv0(ib - m, j, k, n) = 2*svprof(k, n) - sv0(ib - m + 1, j, k, n)
             svm(ib - m, j, k, n) = 2*svprof(k, n) - svm(ib - m + 1, j, k, n)
           end do
         end do
       end do
     end do

   end subroutine xsi_profile


     subroutine xsi_custom
       use modglobal,    only : ib, ie, jb, je, jtot, kb, ke, nsv, ihc, khc
       use modfields,    only : sv0, svm, svprof
       use decomp_2d,    only : zstart

       integer j, k, n, m

       do j = jb, je
          if (j + zstart(2) - 1 == jtot/2) then
             do k = kb, ke + khc
                do n = 1, nsv
                   do m = 1, ihc
                      sv0(ib - m, j-1:j+1, k, n) = 2*svprof(k, n) - sv0(ib - m + 1, j-1:j+1, k, n)
                      svm(ib - m, j-1:j+1, k, n) = 2*svprof(k, n) - svm(ib - m + 1, j-1:j+1, k, n)
                   end do
                end do
             end do
          end if
       end do

   end subroutine xsi_custom


   subroutine xsi_driver
     use modglobal,    only : ib, ie, ihc, jb, je, jhc, kb, ke, khc, nsv
     use modinletdata, only : sv0driver, svmdriver
     use modfields,    only : sv0, svm

     integer j, k, n, m

     do j = jb - jhc, je + jhc
       do k = kb, ke + khc
         do n = 1, nsv
           do m = 1, ihc
             sv0(ib - m, j, k, n) = sv0driver(j, k, n)
             svm(ib - m, j, k, n) = svmdriver(j, k, n)
           end do
         end do
       end do
     end do

   end subroutine xsi_driver


   subroutine xmo_convective
     use modglobal,      only : ie, dxi, rk3step, dt, ih
     use modfields,      only : u0, um, v0, vm, w0, wm, e120, e12m, uouttot
     use modsubgriddata, only : loneeqn
     real rk3coef
     integer i

     rk3coef = dt/(4.-dble(rk3step))

     do i = 1, ih
      v0(ie + i, :, :) = v0(ie+i, :, :) - (v0(ie+i, :, :) - v0(ie+i-1, :, :))*dxi*rk3coef*uouttot
      w0(ie + i, :, :) = w0(ie+i, :, :) - (w0(ie+i, :, :) - w0(ie+i-1, :, :))*dxi*rk3coef*uouttot
      vm(ie + i, :, :) = vm(ie+i, :, :) - (vm(ie+i, :, :) - vm(ie+i-1, :, :))*dxi*rk3coef*uouttot
      wm(ie + i, :, :) = wm(ie+i, :, :) - (wm(ie+i, :, :) - wm(ie+i-1, :, :))*dxi*rk3coef*uouttot

      if (loneeqn) then
        e120(ie + i, :, :) = e120(ie+i-1, :, :) - (e120(ie + i, :, :) - e120(ie+i-1, :, :))*dxi*rk3coef*uouttot
        e12m(ie + i, :, :) = e12m(ie+i-1, :, :) - (e12m(ie + i, :, :) - e12m(ie+i-1, :, :))*dxi*rk3coef*uouttot
      end if
     end do
   end subroutine xmo_convective


   subroutine xmo_Neumann
     use modglobal,      only : ie, ih
     use modfields,      only : u0, um, v0, vm, w0, wm, e120, e12m
     use modsubgriddata, only : loneeqn
     integer i

     do i = 1, ih
      v0(ie + i, :, :) = v0(ie+i-1, :, :)
      w0(ie + i, :, :) = w0(ie+i-1, :, :)
      vm(ie + i, :, :) = vm(ie+i-1, :, :)
      wm(ie + i, :, :) = wm(ie+i-1, :, :)

      if (loneeqn) then
        e120(ie + i, :, :) = e120(ie+i-1, :, :)
        e12m(ie + i, :, :) = e12m(ie+i-1, :, :)
      end if
     end do

   end subroutine xmo_Neumann


   subroutine xTo_convective
     use modglobal, only : ie, dxi, rk3step, dt, ih
     use modfields, only : thl0, thlm, uouttot
     real rk3coef
     integer i

     rk3coef = dt/(4.-dble(rk3step))

     do i = 1, ih
      thl0(ie + i, :, :) = thl0(ie+i, :, :) - (thl0(ie+i, :, :) - thl0(ie+i-1, :, :))*dxi*rk3coef*uouttot
      thlm(ie + i, :, :) = thlm(ie+i, :, :) - (thlm(ie+i, :, :) - thlm(ie+i-1, :, :))*dxi*rk3coef*uouttot
     end do

   end subroutine xTo_convective


   subroutine xTo_Neumann
     use modglobal, only : ie, ih
     use modfields, only : thl0, thlm
     integer i

     do i = 1, ih
      thl0(ie + i, :, :) = thl0(ie+i-1, :, :)
      thlm(ie + i, :, :) = thlm(ie+i-1, :, :)
     end do

   end subroutine xTo_Neumann


   subroutine xqo_convective
     use modglobal, only : ie, dxi, rk3step, dt, ih
     use modfields, only : qt0, qtm, uouttot
     real rk3coef
     integer i

     rk3coef = dt/(4.-dble(rk3step))

     do i = 1, ih
      qt0(ie + i, :, :) = qt0(ie+i-1, :, :) - (qt0(ie+i, :, :) - qt0(ie+i-1, :, :))*dxi*rk3coef*uouttot
      qtm(ie + i, :, :) = qtm(ie+i-1, :, :) - (qtm(ie+i, :, :) - qtm(ie+i-1, :, :))*dxi*rk3coef*uouttot
     end do

   end subroutine xqo_convective


   subroutine xso_convective
     use modglobal, only : ie, rk3step, dt, dxi, nsv, ih
     use modfields, only :sv0, svm, uouttot
     real rk3coef
     integer n, i

     rk3coef = dt/(4.-dble(rk3step))
     do i = 1, ih
     do n = 1, nsv
       sv0(ie + i, :, :, n) = sv0(ie + i, :, :, n) - (sv0(ie + i, :, :, n) - sv0(ie+i-1, :, :, n))*dxi*rk3coef*uouttot
       svm(ie + i, :, :, n) = svm(ie + i, :, :, n) - (svm(ie + i, :, :, n) - svm(ie+i-1, :, :, n))*dxi*rk3coef*uouttot
     end do
     end do

   end subroutine xso_convective


   subroutine xso_Neumann
     use modglobal, only : ie, ihc, rk3step, dt, dxi, nsv
     use modfields, only :sv0, svm
     real rk3coef
     integer n, m

     rk3coef = dt/(4.-dble(rk3step))

     do n = 1, nsv
       do m = 1, ihc
         sv0(ie + m, :, :, n) = sv0(ie, :, :, n)
         svm(ie + m, :, :, n) = svm(ie, :, :, n)
       end do
     end do

   end subroutine xso_Neumann


   subroutine ymi_profile
     use modglobal,      only : ib, ie, jb, je, kb, ke, ih, jh, kh
     use modfields,      only : u0, um, v0, vm, w0, wm, e120, e12m, uprof, vprof, e12prof
     use modsubgriddata, only : loneeqn
     integer i, j, k

     do i = ib - ih, ie + ih
       do k = kb, ke + kh
         v0(i, jb, k) = vprof(k)
         vm(i, jb, k) = vprof(k)
         do j = 1, jh
          v0(i, jb - j, k) = 2*v0(i, jb - j + 1, k) - v0(i, jb - j + 2, k)
          vm(i, jb - j, k) = 2*vm(i, jb - j + 1, k) - vm(i, jb - j + 2, k)
          if (j==1) then
            u0(i, jb - j, k) = 2*uprof(k) - u0(i, jb - j + 1, k)
            um(i, jb - j, k) = 2*uprof(k) - um(i, jb - j + 1, k)
          else
            u0(i, jb - j, k) = 2*u0(i, jb - j + 1, k) - u0(i, jb - j + 2, k)
            um(i, jb - j, k) = 2*um(i, jb - j + 1, k) - um(i, jb - j + 2, k)
          end if
          w0(i, jb - j, k) = -w0(i, jb - j + 1, k)
          wm(i, jb - j, k) = -wm(i, jb - j + 1, k)
         end do
       end do
     end do

     if (loneeqn) then
       do i = ib - ih, ie + ih
         do k = kb, ke + kh
          do j = 1, jh
            if (j==1) then
              e120(i, jb - j, k) = 2*e12prof(k) - e120(i, jb - j + 1, k)
              e12m(i, jb - j, k) = 2*e12prof(k) - e12m(i, jb - j + 1, k)
            else
              e120(i, jb - j, k) = 2*e120(i, jb - j + 1, k) - e120(i, jb - j + 2, k)
              e12m(i, jb - j, k) = 2*e12m(i, jb - j + 1, k) - e12m(i, jb - j + 2, k)
            end if
          end do
         end do
       end do
     end if

   end subroutine ymi_profile


   subroutine yTi_profile
     use modglobal, only : ib, ie, jb, je, kb, ke, ih, jh, kh
     use modfields, only : thl0, thlm, thlprof

     integer i, j, k

     do i = ib - ih, ie + ih
       do k = kb, ke + kh
        do j = 1, jh
          if (j==1) then
            thl0(i, jb - j, k) = 2*thlprof(k) - thl0(i, jb - j + 1, k)
            thlm(i, jb - j, k) = 2*thlprof(k) - thlm(i, jb - j + 1, k)
          else
            thl0(i, jb - j, k) = 2*thl0(i, jb - j + 1, k) - thl0(i, jb - j + 2, k)
            thlm(i, jb - j, k) = 2*thlm(i, jb - j + 1, k) - thlm(i, jb - j + 2, k)
          end if
        end do
       end do
     end do

   end subroutine yTi_profile


   subroutine yqi_profile
     use modglobal, only : ib, ie, jb, je, kb, ke, ih, jh, kh
     use modfields, only : qt0, qtm, qtprof

     integer i, j, k

     do i = jb - jh, ie + jh
       do k = kb, ke + kh
        do j = 1, jh
          if (j==1) then
            qt0(i, jb - j, k) = 2*qtprof(k) - qt0(i, jb - j + 1, k)
            qtm(i, jb - j, k) = 2*qtprof(k) - qtm(i, jb - j + 1, k)
          else
            qt0(i, jb - j, k) = 2*qt0(i, jb - j + 1, k) - qt0(i, jb - j + 2, k)
            qtm(i, jb - j, k) = 2*qtm(i, jb - j + 1, k) - qtm(i, jb - j + 2, k)
          end if
        end do
       end do
     end do

   end subroutine yqi_profile


   subroutine ysi_profile
     use modglobal, only : ib, ie, jb, je, kb, ke, nsv, ihc, jhc, khc
     use modfields, only : sv0, svm, svprof

     integer i, k, n, m

     do i = ib - ihc, ie + ihc
       do k = kb, ke + khc
         do n = 1, nsv
           do m = 1, jhc
             sv0(i, jb - m, k, n) = 2*svprof(k, n) - sv0(i, jb - m + 1, k, n)
             svm(i, jb - m, k, n) = 2*svprof(k, n) - svm(i, jb - m + 1, k, n)
           end do
         end do
       end do
     end do

   end subroutine ysi_profile


   subroutine ymo_convective
     use modglobal,      only : je, dyi, rk3step, dt, jh
     use modfields,      only : u0, um, v0, vm, w0, wm, e120, e12m, vouttot
     use modsubgriddata, only : loneeqn

     real rk3coef
     integer j

     rk3coef = dt/(4.-dble(rk3step))

     do j = 1, jh
      ! change to vouttot
      u0(:, je + j, :) = u0(:, je + j, :) - (u0(:, je + j, :) - u0(:, je + j - 1, :))*dyi*rk3coef*vouttot
      um(:, je + j, :) = um(:, je + j, :) - (um(:, je + j, :) - um(:, je + j - 1, :))*dyi*rk3coef*vouttot
      w0(:, je + j, :) = w0(:, je + j, :) - (w0(:, je + j, :) - w0(:, je + j - 1, :))*dyi*rk3coef*vouttot
      wm(:, je + j, :) = wm(:, je + j, :) - (wm(:, je + j, :) - wm(:, je + j - 1, :))*dyi*rk3coef*vouttot

      if (loneeqn) then
        e120(:, je + j, :) = e120(:, je + j, :) - (e120(:, je + j, :) - e120(:, je + j - 1, :))*dyi*rk3coef*vouttot
        e12m(:, je + j, :) = e12m(:, je + j, :) - (e12m(:, je + j, :) - e12m(:, je + j - 1, :))*dyi*rk3coef*vouttot
      end if
     end do

   end subroutine ymo_convective


   subroutine yTo_convective

     use modglobal, only : je, dyi, rk3step, dt, jh
     use modfields, only : thl0, thlm, v0, vouttot

     real rk3coef
     integer j

     rk3coef = dt/(4.-dble(rk3step))

     do j = 1, jh
      thl0(:, je + j, :) = thl0(:, je + j, :) - (thl0(:, je + j, :) - thl0(:, je + j - 1, :))*dyi*rk3coef*vouttot
      thlm(:, je + j, :) = thlm(:, je + j, :) - (thlm(:, je + j, :) - thlm(:, je + j - 1, :))*dyi*rk3coef*vouttot
     end do

   end subroutine yTo_convective


   subroutine yqo_convective

     use modglobal, only : je, dyi, rk3step, dt, jh
     use modfields, only : qt0, qtm, v0, vouttot

     real rk3coef
     integer j
     rk3coef = dt/(4.-dble(rk3step))

     do j = 1, jh
      qt0(:, je + j, :) = qt0(:, je + j, :) - (qt0(:, je + j, :) - qt0(:, je + j - 1, :))*dyi*rk3coef*vouttot
      qtm(:, je + j, :) = qtm(:, je + j, :) - (qtm(:, je + j, :) - qtm(:, je + j - 1, :))*dyi*rk3coef*vouttot
     end do

   end subroutine yqo_convective


   subroutine yso_convective

     use modglobal, only : je, rk3step, dt, dyi, nsv, jhc
     use modfields, only :sv0, svm, v0, vouttot

     real rk3coef
     integer n, j

     rk3coef = dt/(4.-dble(rk3step))

     do j = 1, jhc
     do n = 1, nsv
       sv0(:, je + j, :, n) = sv0(:, je + j, :, n) - (sv0(:, je + j, :, n) - sv0(:, je + j - 1, :, n))*dyi*rk3coef*vouttot
       svm(:, je + j, :, n) = svm(:, je + j, :, n) - (svm(:, je + j, :, n) - svm(:, je + j - 1, :, n))*dyi*rk3coef*vouttot
     end do
     end do

   end subroutine yso_convective


   subroutine yso_Neumann

     use modglobal, only : je, jhc, rk3step, dt, dyi, nsv
     use modfields, only : sv0, svm

     real rk3coef
     integer n, m

     rk3coef = dt/(4.-dble(rk3step))

     do n = 1, nsv
       do m = 1, jhc
         sv0(:, je + m, :, n) = sv0(:, je, :, n)
         svm(:, je + m, :, n) = svm(:, je, :, n)
       end do
     end do

   end subroutine yso_Neumann


   !>set boundary conditions pup,pvp,pwp in subroutine fillps in modpois.f90
   subroutine bcpup(pup, pvp, pwp, rk3coef)

     use modglobal,    only : ib, ie, jb, je, ih, jh, kb, ke, kh, rk3step, dxi, dyi, dzhi, &
                              ibrank, ierank, jbrank, jerank, BCxm, BCym, BCtopm, &
                              BCtopm_freeslip, BCtopm_noslip, BCtopm_pressure, &
                              BCxm_periodic, BCxm_profile, BCxm_driver, &
                              BCym_periodic, BCym_profile
     use modfields,    only : pres0, up, vp, wp, um, vm, wm, w0, u0, v0, uouttot, vouttot, uinit, vinit, uprof, vprof, pres0, IIc, IIcs
     use modmpi,       only : excjs, excis, myid, avexy_ibm
     use modinletdata, only : u0driver
     use decomp_2d,    only : exchange_halo_z

     real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb:ke + kh), intent(inout) :: pup
     real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb:ke + kh), intent(inout) :: pvp
     real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb:ke + kh), intent(inout) :: pwp
     real, dimension(kb:ke+kh) :: pres0ij

     real, intent(in) :: rk3coef
     real rk3coefi

     integer i, j, k

     rk3coefi = 1./rk3coef

     ! if (jbrank) write(*,*) "jb before exhange_halo ", pvp(ie/2,jb,ke)
     ! if (jerank) write(*,*) "je before exhange_halo ", pvp(ie/2,je+1,ke)
     ! Watch this communication as it is slightly different to normal -
     ! maybe safer to just resize to kb-kh:ke+kh
     call exchange_halo_z(pup, opt_zlevel=(/ih,jh,0/))
     call exchange_halo_z(pvp, opt_zlevel=(/ih,jh,0/))
     call exchange_halo_z(pwp, opt_zlevel=(/ih,jh,0/))
     ! if (jbrank) write(*,*) "jb after exhange_halo ", pvp(ie/2,jb,ke)
     ! if (jerank) write(*,*) "je after exhange_halo ", pvp(ie/2,je+1,ke)

     select case(BCtopm)
     case(BCtopm_freeslip, BCtopm_noslip)
       do j = jb, je
         do i = ib, ie
           pwp(i, j, kb) = 0.
           pwp(i, j, ke + kh) = 0.
         end do
       end do

     case(BCtopm_pressure)
       call avexy_ibm(pres0ij(kb:ke+kh),pres0(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.false.)

       do j = jb, je
         do i = ib, ie
           pwp(i, j, kb)  = 0.
           !pwp(i, j, ke + 1) = wm(i, j, ke+1) * rk3coefi - (-pres0ij(ke) - pres0(i,j,ke)) * dzhi(ke+1) ! Doesn't work
           do k = 1, kh
            pwp(i, j, ke + k) = wm(i, j, ke+k) * rk3coefi + 2 * pres0ij(ke+k-1)*dzhi(ke+k)
            wp(i, j, ke + k) = pwp(i, j, ke+k) - wm(i,j,ke+k) * rk3coefi
           end do
         end do
       end do
     end select !BCtopm

     select case(BCxm)
     case(BCxm_periodic)
       if (ibrank .and. ierank) then ! not parallelised in x
         do k = kb, ke
            do j = jb, je
              do i = 1, ih
               pup(ie+i, j, k) = pup(ib-i+1, j, k) ! cyclic
              end do
            end do
         end do
       end if

     case(BCxm_profile)
       if (ibrank) then
         do k=kb,ke
           do j=jb-jh,je+jh
             pup(ib, j, k) = uprof(k) * rk3coefi
             up(ib, j, k) = 0.
           end do
         end do
       end if

       if (ierank) then
         do k = kb+kh, ke
           do j = jb-jh, je+jh
            do i = 1, ih
             ! convective
             pup(ie+i, j, k) = um(ie+i, j, k) * rk3coefi - (u0(ie+i, j, k) - u0(ie+i-1,j,k)) * dxi * uouttot !u0(ie,j,k) ! du/dt +u*du/dx=0 -> pup(i)=um(i)/rk3coef -um(i)*(um(i)-um(i-1))/dxf(i-1)
             ! Neumann
             !pup(ie+1,j,k) = pup(ie,j,k)
             up(ie+i, j, k) = pup(ie+i, j, k) - um(ie+i,j,k)*rk3coefi
            end do
           end do
         end do

         do i = 1, ih
          ! Neumann at bottom - performs better with no slip
          pup(ie+i, :, kb) = pup(ie+i-1, :, kb)
          up(ie+i, :, kb) = pup(ie+i,: , kb) - um(ie+i, :, kb) * rk3coefi
         end do
       end if

     case(BCxm_driver)
       if (ibrank) then
         do k = kb, ke
           do j = jb - jh, je + jh
             pup(ib, j, k) = u0driver(j, k) * rk3coefi
             up(ib, j, k) = 0. ! u(ib) only evolves according to pressure correction
           end do
         end do
       end if

       if (ierank) then
         do k = kb, ke
           do j = jb-jh, je+jh
            do i = 1, ih
             pup(ie+i, j, k) = um(ie+i, j, k) * rk3coefi - (u0(ie+i, j, k) - u0(ie+i-1, j, k)) * dxi * uouttot    ! du/dt +u*du/dx=0 -> pup(i)=um(i)/rk3coef -um(i)*(um(i)-um(i-1))/dxf(i-1)
             ! !Neumann
             !pup(ie+1,j,k) = pup(ie,j,k)
             up(ie+i, j, k) = pup(ie+i, j, k) - um(ie+i, j, k) * rk3coefi
            end do
           end do
         end do
         ! Neumann at bottom - performs better with no slip
         ! pup(ie+1, :, kb) = pup(ie, :, kb)
         ! up(ie+1, :, kb) = pup(ie+1, :, kb) - um(ie+1, :, kb) * rk3coefi
       end if
    end select ! BCxm

    select case(BCym)
    case(BCym_periodic)
      if (jbrank .and. jerank) then ! not parallelised in y
        do k = kb, ke
           do i = ib, ie
            do j = 1, jh
              pvp(i, je+j, k) = pvp(i, jb-j+1, k) ! cyclic
            end do
           end do
        end do
      end if

    case(BCym_profile)
      if (jbrank) then
        do k = kb, ke
          do i = ib-ih, ie+ih
            pvp(i,jb,k) = vprof(k)*rk3coefi
            vp(i,jb,k) = 0.
          end do
        end do
      end if

      if (jerank) then
        do k = kb, ke
          do i = ib-ih, ie+ih
            do j = 1, jh
              ! change to vouttot
              pvp(i, je+j, k) = vm(i, je+j, k) * rk3coefi - (v0(i, je+j, k) - v0(i, je+j-1, k)) * dyi * vouttot
              vp(i, je+j, k) = pvp(i, je+j, k) - vm(i,je+j,k)*rk3coefi
            end do
          end do
        end do
        do j = 1, jh
          pvp(:, je+j, kb) = pvp(:, je+j-1, kb)
          vp(:, je+j, kb) = pvp(:, je+j, kb) - vm(:, je+j, kb)*rk3coefi
        end do
      end if

    end select

   end subroutine bcpup

   !>set pressure boundary conditions
   subroutine bcp(p)

     use modglobal, only : ib, ie, jb, je, ih, jh, kb, ke, kh, dyi, rk3step, dt, &
                           ibrank, ierank, jbrank, jerank, BCxm, BCym, BCxm_periodic, BCym_periodic
     use modfields, only : pres0, up, u0, um, uouttot, vp, v0
     use decomp_2d, only : exchange_halo_z

     real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb - kh:ke + kh), intent(inout) :: p !< pressure
     integer i, j, k
     real rk3coef, rk3coefi

     if (rk3step == 0) then ! dt not defined yet
       rk3coef = 1.
     else
       rk3coef = dt / (4. - dble(rk3step))
     end if
     rk3coefi = 1. / rk3coef

     call exchange_halo_z(p)
     call exchange_halo_z(pres0)

     if (BCxm .eq. BCxm_periodic) then
       if (ibrank .and. ierank) then
         do j = jb, je
           do k = kb, ke
            do i = 1, ih
             p(ib-i, j, k) = p(ie+i-1, j, k)
             p(ie+i, j, k) = p(ib-i+1, j, k)
             !pres0(ib - 1, j, k) = pres0(ie, j, k)
             !pres0(ie + 1, j, k) = pres0(ib, j, k)
            end do
           end do
         end do
       end if

     else
       if (ibrank) then
         do k = kb, ke
           do j = jb-jh, je+jh
            do i = 1, ih
             p(ib-i, j, k) = p(ib-i+1, j, k)
             pres0(ib-i, j, k) = pres0(ib-i+1, j, k)
            end do
           end do
         end do
       end if

       if (ierank) then
         do k = kb, ke
           do j = jb-jh, je+jh
            do i = 1, ih
             p(ie+i, j, k) = p(ie+i-1, j, k)
             pres0(ie+i, j, k) = pres0(ie+i-1, j, k)
            end do
           end do
         end do
       end if

     end if ! BCxm

     if (BCym .eq. BCym_periodic) then
       if (jbrank .and. jerank) then
         do i = ib, ie
           do k = kb, ke
            do j = 1, jh
             p(i, jb-j, k) = p(i, je+j-1, k)
             p(i, je+j, k) = p(i, jb-j+1, k)
             !pres0(ib - 1, j, k) = pres0(ie, j, k)
             !pres0(ie + 1, j, k) = pres0(ib, j, k)
            end do
           end do
         end do
       end if
     else
       if (jbrank) then
         do k = kb, ke
           do i = ib-ih, ie+ih
            do j = 1, jh
             p(i,jb-j,k) = p(i,jb-j+1,k)
             pres0(i,jb-j,k) = pres0(i,jb-j+1,k)
            end do
           enddo
         enddo
       end if

       if (jerank) then
         do k = kb, ke
           do i = ib-ih, ie+ih
            do j = 1, jh
             p(i, je+j, k) = p(i,je+j-1,k)
             pres0(i, je+j, k) = pres0(i,je+j-1,k)
            end do
           end do
         end do
       end if

     end if !BCym

   end subroutine bcp

   !>
   !! grwdamp damps gravity waves in the upper part of the domain.
   !>
   !! The lower limit of the damping region is set by ksp
   !! Horizontal fluctuations at the top of the domain (for instance gravity waves)
   !! are damped out by a sponge layer through an additional forcing/source term.
   !! \latexonly
   !! \begin{eqnarray}
   !! \force{i}{sp}(z) &=& -\frac{1}{t^{\mr{sp}}}\left(\xav{\fav{u_i}}-\fav{u_i}\right), \\\\
   !!  \source{\varphi}{sp}(z) &=& -\frac{1}{t^{\mr{sp}}}\left(\xav{\varphi}-\varphi\right),
   !! \end{eqnarray}
   !! with $t^{\mr{sp}}$ a relaxation time scale that goes from
   !! $t^{\mr{sp}}_0=1/(2.75\times10^{-3})\mr{s}\approx 6$min at the top of the domain
   !! to infinity at the bottom of the sponge layer.
   !! \endlatexonly
   subroutine grwdamp
      use modglobal, only:ke, kmax, lcoriol, igrw_damp, geodamptime
      use modfields, only:up, vp, wp, thlp, qtp, u0, v0, w0, thl0, qt0, ug, vg, thl0av, qt0av, u0av, v0av
      use modmpi, only:myid
      implicit none

      integer k
      select case (igrw_damp)
      case (0) !do nothing
      case (1)
         do k = ksp, ke
            up(:, :, k) = up(:, :, k) - (u0(:, :, k) - u0av(k))*tsc(k)
            vp(:, :, k) = vp(:, :, k) - (v0(:, :, k) - v0av(k))*tsc(k)
            wp(:, :, k) = wp(:, :, k) - w0(:, :, k)*tsc(k)
            thlp(:, :, k) = thlp(:, :, k) - (thl0(:, :, k) - thl0av(k))*tsc(k)
            qtp(:, :, k) = qtp(:, :, k) - (qt0(:, :, k) - qt0av(k))*tsc(k)
         end do
         if (lcoriol) then
            do k = ksp, ke
               up(:, :, k) = up(:, :, k) - (u0(:, :, k) - ug(k))*((1./(geodamptime*rnu0))*tsc(k))
               vp(:, :, k) = vp(:, :, k) - (v0(:, :, k) - vg(k))*((1./(geodamptime*rnu0))*tsc(k))
            end do
         end if
      case (2)
         do k = ksp, ke
            up(:, :, k) = up(:, :, k) - (u0(:, :, k) - ug(k))*tsc(k)
            vp(:, :, k) = vp(:, :, k) - (v0(:, :, k) - vg(k))*tsc(k)
            wp(:, :, k) = wp(:, :, k) - w0(:, :, k)*tsc(k)
            thlp(:, :, k) = thlp(:, :, k) - (thl0(:, :, k) - thl0av(k))*tsc(k)
            qtp(:, :, k) = qtp(:, :, k) - (qt0(:, :, k) - qt0av(k))*tsc(k)
         end do
      case (3)
         do k = ksp, ke
            up(:, :, k) = up(:, :, k) - (u0(:, :, k) - u0av(k))*tsc(k)
            vp(:, :, k) = vp(:, :, k) - (v0(:, :, k) - v0av(k))*tsc(k)
            wp(:, :, k) = wp(:, :, k) - w0(:, :, k)*tsc(k)
            thlp(:, :, k) = thlp(:, :, k) - (thl0(:, :, k) - thl0av(k))*tsc(k)
            qtp(:, :, k) = qtp(:, :, k) - (qt0(:, :, k) - qt0av(k))*tsc(k)
         end do
      case default
         write(0, *) "ERROR: no gravity wave damping option selected"
         stop 1
      end select

      return
   end subroutine grwdamp


   subroutine fluxtop(field, ek, flux)
      use modglobal, only:ib, ie, ih, jb, je, jh, kb, ke, kh, dzf, dzh, dzhi, eps1

      real, intent(inout) :: field(ib - ih:ie + ih, jb - jh:je + jh, kb - kh:ke + kh)
      real, intent(in)    ::    ek(ib - ih:ie + ih, jb - jh:je + jh, kb - kh:ke + kh)
      real, intent(in)    :: flux
      integer k
      !
      if (abs(flux) .le. eps1) then !it's zero-flux, we don't need to do the calculation
         do k = 1, kh
          field(:, :, ke + k) = field(:, :, ke + k - 1)
         end do
      else
        do k = 1, kh
         field(:, :, ke + k) = field(:, :, ke + k - 1) + dzh(ke + k)*flux/(dzhi(ke + k)*(0.5*(dzf(ke+k-1)*ek(:, :, ke + k) + dzf(ke + k)*ek(:, :, ke+k-1))))
        end do
      end if
      !
   end subroutine fluxtop

   subroutine valuetop(field, val)
      use modglobal, only:ib, ie, ih, jb, je, jh, kb, ke, kh, dzh, dzf, dzhi, dzfi
      use modmpi, only : myid
      real, intent(inout) :: field(ib - ih:ie + ih, jb - jh:je + jh, kb - kh:ke + kh)
      real, intent(in)    :: val
      integer k

      ! (field(i, j, kp)*dzf(k) + field(i, j, k)*dzf(kp))*dzhi(kp)*0.5 = val
      !field(:,:,ke+1) = (2.*val*dzh(ke+1) - field(:,:,ke)*dzf(ke+1)) * dzfi(ke)
      do k = 1, kh
        field(:, :, ke + k) = 2*val - field(:, :, ke + k - 1)
      end do
      !if (myid == 0) write(*,*) (field(40, 1, ke+1)*dzf(ke) + field(40, 1, ke)*dzf(ke+1))*dzhi(ke+1)*0.5

   end subroutine valuetop

   subroutine fluxtopscal(flux)
      use modglobal, only:ib, ie, ih, jb, je, jh, kb, ke, kh, dzf, dzh, dzhi, nsv, khc
      use modfields, only:sv0, svm
      use modsubgriddata, only:ekh

      real, intent(in)    :: flux(1:nsv)
      integer :: m, n
      !
      !all the ghost cells have the same value?
      do m = 1, khc
      do n = 1, nsv
  sv0(ib-ih:ie+ih,jb-jh:je+jh,ke+m,n) = sv0(ib-ih:ie+ih,jb-jh:je+jh,ke,n) + dzh(ke+1) * flux(n) / ( dzhi(ke+1) * (0.5*(dzf(ke)*ekh(ib-ih:ie+ih,jb-jh:je+jh,ke+1)+dzf(ke+1)*ekh(ib-ih:ie+ih,jb-jh:je+jh,ke))))
  svm(ib-ih:ie+ih,jb-jh:je+jh,ke+m,n) = svm(ib-ih:ie+ih,jb-jh:je+jh,ke,n) + dzh(ke+1) * flux(n) / ( dzhi(ke+1) * (0.5*(dzf(ke)*ekh(ib-ih:ie+ih,jb-jh:je+jh,ke+1)+dzf(ke+1)*ekh(ib-ih:ie+ih,jb-jh:je+jh,ke))))
      end do
      end do
      !
   end subroutine fluxtopscal

   subroutine valuetopscal(val)
      use modglobal, only:ib, ie, ih, jb, je, jh, kb, ke, kh, eps1, nsv, khc
      use modfields, only:sv0, svm
      real, intent(in)    :: val(1:nsv)
      integer :: m, n
      !
      ! all the ghost cells have the same vlaue?
      do m = 1, khc
      do n = 1, nsv
         sv0(: , : , ke + m, n) = 2*val(n) - sv0(: , : , ke, n)
         svm(: , : , ke + m, n) = 2*val(n) - svm(: , : , ke, n)
      end do
      end do
      !
   end subroutine valuetopscal


!>Set thl, qt and sv(n) equal to slab average at level kmax
! Think this can be removed
   subroutine tqaver

      use modmpi, only:comm3d, mpierr, my_real, mpi_sum
      use modglobal, only:ib, ie, jb, je, ih, jh, kb, ke, nsv, rslabs
      use modfields, only:thl0, qt0, sv0
      implicit none

      real thl0a, qt0a
      real thl0al, qt0al
      integer n
      real, allocatable, dimension(:) :: sv0al, sv0a
      allocate (sv0al(nsv), sv0a(nsv))

      thl0al = sum(thl0(ib:ie, jb:je, ke))
      qt0al = sum(qt0(ib:ie, jb:je, ke))

      do n = 1, nsv
         sv0al(n) = sum(sv0(ib:ie, jb:je, ke, n))
      enddo

      call MPI_ALLREDUCE(thl0al, thl0a, 1, MY_REAL, &
                         MPI_SUM, comm3d, mpierr)
      call MPI_ALLREDUCE(qt0al, qt0a, 1, MY_REAL, &
                         MPI_SUM, comm3d, mpierr)
      if (nsv > 0) then
         call MPI_ALLREDUCE(sv0al, sv0a, nsv, MY_REAL, &
                            MPI_SUM, comm3d, mpierr)
      end if

      thl0a = thl0a/rslabs
      qt0a = qt0a/rslabs
      sv0a = sv0a/rslabs

      thl0(ib:ie, jb:je, ke) = thl0a
      qt0(ib:ie, jb:je, ke) = qt0a
      do n = 1, nsv
         sv0(ib:ie, jb:je, ke, n) = sv0a(n)
      enddo
      deallocate (sv0al, sv0a)

      return
   end subroutine tqaver

end module modboundary
