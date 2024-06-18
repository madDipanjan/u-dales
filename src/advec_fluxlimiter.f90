!> \file advec_fluxlimiter.f90
!!  Does advection with a flux limiter scheme between 2nd order central differencing  and 1st order upwind.
!! \par Revision list
!! \par Authors
!! Second order central differencing can be used for variables where neither very
!! high accuracy nor strict monotonicity is necessary.
!! \latexonly
!!\begin{eqnarray}
!! F_{i-\frac{1}{2}}^{2nd} &=&
!!\fav{u}_{i-\frac{1}{2}}\frac{\phi_{i}+\phi_{i-1}}{2},
!!\end{eqnarray}
!! \endlatexonly
!!
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

!> Advection at cell center
subroutine advecc_Flimiter(hi, hj, hk, putin, putout)
    use modglobal, only:ib, ie, jb, je, kb, ke, dxi, dyi, dzfi, dzh, dzf
    use modfields, only:u0, v0, w0
    implicit none
 
    integer, intent(in) :: hi !< size of halo in i
    integer, intent(in) :: hj !< size of halo in j
    integer, intent(in) :: hk !< size of halo in k
    real, dimension(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk), intent(in)  :: putin !< Input: the cell centered field
    real, dimension(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk), intent(inout) :: putout !< Output: the tendency
 
    real :: flux_xL, flux_xR, flux_yF, flux_yB, flux_zT, flux_zB
    integer :: i, j, k, ip, im, jp, jm, kp, km, ip2, im2, jp2, jm2, kp2, km2

    do k = kb, ke
        km = k - 1
        kp = k + 1
        km2 = k - 2
        kp2 = k + 2
        do j = jb, je
            jm = j - 1
            jp = j + 1
            jm2 = j - 2
            jp2 = j + 2
            do i = ib, ie
             im = i - 1
             ip = i + 1
             im2 = i - 2
             ip2 = i + 2

             call computefluxlimiter(putin(ip2,j,k), putin(ip,j,k), putin(i,j,k), putin(im,j,k), putin(im2,j,k), &
                                     u0(i,j,k), u0(ip,j,k), flux_xL,flux_xR)

             call computefluxlimiter(putin(i,jp2,k), putin(i,jp,k), putin(i,j,k), putin(i,jm,k), putin(i,jm2,k), &
                                     v0(i,j,k), v0(i,jp,k), flux_yF,flux_yB)

             if (k==kb) then
                call compFLimiter_nonuni(putin(i,j,kp2), putin(i,j,kp), putin(i,j,k), putin(i,j,km), putin(i,j,km2), &
                                         dzh(kp2), dzh(kp), dzh(k), dzh(k), dzf(kp), dzf(k), dzf(km), &
                                         w0(i,j,k), w0(i,j,kp), flux_zB,flux_zT)
             else
                call compFLimiter_nonuni(putin(i,j,kp2), putin(i,j,kp), putin(i,j,k), putin(i,j,km), putin(i,j,km2), &
                                         dzh(kp2), dzh(kp), dzh(k), dzh(km), dzf(kp), dzf(k), dzf(km), &
                                         w0(i,j,k), w0(i,j,kp), flux_zB,flux_zT)
             end if

             putout(i, j, k) = putout(i, j, k) - ( &
                               ( flux_xR - flux_xL )*dxi &                ! d(uc)/dx
                             + ( flux_yB - flux_yF )*dyi &                ! d(vc)/dy
                             + ( flux_zT - flux_zB )*dzfi(k) )            ! d(wc)/dz
            end do
        end do
    end do
end subroutine advecc_Flimiter


!> Advection at the u point.
subroutine advecu_Flimiter(putin, putout)
    use modglobal, only:ib, ie, ih, jb, je, jh, kb, ke, kh, dxi, dxi5, dyi5, dzfi5, dzh, dzf
    use modfields, only:u0, v0, w0, pres0
    implicit none

    real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb - kh:ke + kh), intent(in) :: putin !< Input: the u-field
    real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb:ke + kh), intent(inout) :: putout !< Output: the tendency

    real :: flux_xL, flux_xR, flux_yF, flux_yB, flux_zT, flux_zB
    integer :: i, j, k, ip, im, jp, jm, kp, km, ip2, im2, jp2, jm2, kp2, km2

    do k = kb, ke
        km = k - 1
        kp = k + 1
        km2 = k - 2
        kp2 = k + 2
        do j = jb, je
            jm = j - 1
            jp = j + 1
            jm2 = j - 2
            jp2 = j + 2
            do i = ib, ie
              im = i - 1
              ip = i + 1
              im2 = i - 2
              ip2 = i + 2

              call computefluxlimiter(putin(ip2,j,k), putin(ip,j,k), putin(i,j,k), putin(im,j,k), putin(im2,j,k), &
                                      (u0(i,j,k)+u0(im,j,k)), (u0(ip,j,k)+u0(i,j,k)), flux_xL,flux_xR)
          
              call computefluxlimiter(putin(i,jp2,k), putin(i,jp,k), putin(i,j,k), putin(i,jm,k), putin(i,jm2,k), &
                                      (v0(i,j,k)+v0(im,j,k)), (v0(i,jp,k)+v0(im,jp,k)), flux_yF,flux_yB)

              if (k==kb) then
                call compFLimiter_nonuni(putin(i,j,kp2), putin(i,j,kp), putin(i,j,k), putin(i,j,km), putin(i,j,km2), &
                                         dzh(kp2), dzh(kp), dzh(k), dzh(k), dzf(kp), dzf(k), dzf(km), &
                                         (w0(i,j,k)+w0(im,j,k)), (w0(i,j,kp)+w0(im,j,kp)), flux_zB,flux_zT)
              else
                call compFLimiter_nonuni(putin(i,j,kp2), putin(i,j,kp), putin(i,j,k), putin(i,j,km), putin(i,j,km2), &
                                         dzh(kp2), dzh(kp), dzh(k), dzh(km), dzf(kp), dzf(k), dzf(km), &
                                         (w0(i,j,k)+w0(im,j,k)), (w0(i,j,kp)+w0(im,j,kp)), flux_zB,flux_zT)
              end if

              putout(i, j, k) = putout(i, j, k) - ( &
                                ( flux_xR - flux_xL )*dxi5 &                ! d(uu)/dx
                              + ( flux_yB - flux_yF )*dyi5 &                ! d(vu)/dy
                              + ( flux_zT - flux_zB )*dzfi5(k) ) &          ! d(wu)/dz
                              - ((pres0(i, j, k) - pres0(im, j, k))*dxi)    ! - dp/dx
            end do
        end do
    end do
end subroutine advecu_Flimiter


!> Advection at the v point.
subroutine advecv_Flimiter(putin, putout)
    use modglobal, only:ib, ie, ih, jb, je, jh, kb, ke, kh, dyi, dxi5, dyi5, dzfi5, dzh, dzf
    use modfields, only:u0, v0, w0, pres0
    implicit none

    real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb - kh:ke + kh), intent(in) :: putin !< Input: the u-field
    real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb:ke + kh), intent(inout) :: putout !< Output: the tendency
 
    real :: flux_xL, flux_xR, flux_yF, flux_yB, flux_zT, flux_zB
    integer :: i, j, k, ip, im, jp, jm, kp, km, ip2, im2, jp2, jm2, kp2, km2
    
    do k = kb, ke
        km = k - 1
        kp = k + 1
        km2 = k - 2
        kp2 = k + 2
        do j = jb, je
            jm = j - 1
            jp = j + 1
            jm2 = j - 2
            jp2 = j + 2
            do i = ib, ie
             im = i - 1
             ip = i + 1
             im2 = i - 2
             ip2 = i + 2

             call computefluxlimiter(putin(ip2,j,k), putin(ip,j,k), putin(i,j,k), putin(im,j,k), putin(im2,j,k), &
                                     (u0(i,j,k)+u0(i,jm,k)), (u0(ip,j,k)+u0(ip,jm,k)), flux_xL,flux_xR)
                                    
             call computefluxlimiter(putin(i,jp2,k), putin(i,jp,k), putin(i,j,k), putin(i,jm,k), putin(i,jm2,k), &
                                     (v0(i,jm,k)+v0(i,j,k)), (v0(i,j,k)+v0(i,jp,k)), flux_yF,flux_yB)

             if (k==kb) then
                call compFLimiter_nonuni(putin(i,j,kp2), putin(i,j,kp), putin(i,j,k), putin(i,j,km), putin(i,j,km2), &
                                         dzh(kp2), dzh(kp), dzh(k), dzh(k), dzf(kp), dzf(k), dzf(km), &
                                         (w0(i,j,k)+w0(i,jm,k)), (w0(i,j,kp)+w0(i,jm,kp)), flux_zB,flux_zT)
             else
                call compFLimiter_nonuni(putin(i,j,kp2), putin(i,j,kp), putin(i,j,k), putin(i,j,km), putin(i,j,km2), &
                                         dzh(kp2), dzh(kp), dzh(k), dzh(km), dzf(kp), dzf(k), dzf(km), &
                                         (w0(i,j,k)+w0(i,jm,k)), (w0(i,j,kp)+w0(i,jm,kp)), flux_zB,flux_zT)
             end if
                    
             putout(i, j, k) = putout(i, j, k) - ( &
                               ( flux_xR - flux_xL )*dxi5 &                ! d(uv)/dx
                             + ( flux_yB - flux_yF )*dyi5 &                ! d(vv)/dy
                             + ( flux_zT - flux_zB )*dzfi5(k) ) &          ! d(wv)/dz
                             - ((pres0(i, j, k) - pres0(i, jm, k))*dyi)    ! - dp/dy
            end do
        end do
    end do
end subroutine advecv_Flimiter


 !> Advection at the w point.
subroutine advecw_Flimiter(putin, putout)
    use modglobal, only:ib, ie, ih, jb, je, jh, kb, ke, kh, dzhi, dxi5, dyi5, dzh, dzf
    use modfields, only:u0, v0, w0, pres0
    implicit none

    real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb - kh:ke + kh), intent(in) :: putin !< Input: the u-field
    real, dimension(ib - ih:ie + ih, jb - jh:je + jh, kb:ke + kh), intent(inout) :: putout !< Output: the tendency
 
    real :: flux_xL, flux_xR, flux_yF, flux_yB, flux_zT, flux_zB
    integer :: i, j, k, ip, im, jp, jm, kp, km, ip2, im2, jp2, jm2, kp2, km2
 
    do k = kb + 1, ke
        km = k - 1
        kp = k + 1
        km2 = k - 2
        kp2 = k + 2
        do j = jb, je
            jm = j - 1
            jp = j + 1
            jm2 = j - 2
            jp2 = j + 2
            do i = ib, ie
             im = i - 1
             ip = i + 1
             im2 = i - 2
             ip2 = i + 2
 
             call computefluxlimiter(putin(ip2,j,k), putin(ip,j,k), putin(i,j,k), putin(im,j,k), putin(im2,j,k), &
                                     (dzf(km)*u0(i,j,k) + dzf(k)*u0(i,j,km))*dzhi(k), (dzf(km)*u0(ip,j,k) + dzf(k)*u0(ip,j,km))*dzhi(k), &
                                     flux_xL,flux_xR)
                                    
             call computefluxlimiter(putin(i,jp2,k), putin(i,jp,k), putin(i,j,k), putin(i,jm,k), putin(i,jm2,k), &
                                     (dzf(km)*v0(i,j,k) + dzf(k)*v0(i,j,km))*dzhi(k), (dzf(km)*v0(i,jp,k) + dzf(k)*v0(i,jp,km))*dzhi(k), &
                                     flux_yF,flux_yB)

             if (k==kb) then
                call compFLimiter_nonuni(putin(i,j,kp2), putin(i,j,kp), putin(i,j,k), putin(i,j,km), putin(i,j,km2), &
                                         dzf(kp), dzf(k), dzf(km), dzf(km2), dzh(kp), dzh(k), dzh(k), &
                                         (w0(i,j,km)+w0(i,j,k)), (w0(i,j,k)+w0(i,j,kp)), flux_zB,flux_zT)
             else
                call compFLimiter_nonuni(putin(i,j,kp2), putin(i,j,kp), putin(i,j,k), putin(i,j,km), putin(i,j,km2), &
                                         dzf(kp), dzf(k), dzf(km), dzf(km2), dzh(kp), dzh(k), dzh(km), &
                                         (w0(i,j,km)+w0(i,j,k)), (w0(i,j,k)+w0(i,j,kp)), flux_zB,flux_zT)
             end if

             putout(i, j, k) = putout(i, j, k) - ( &
                               ( flux_xR - flux_xL )*dxi5 &                     ! d(uw)/dx
                             + ( flux_yB - flux_yF )*dyi5 &                     ! d(vw)/dy
                             + ( flux_zT - flux_zB )*dzhi(k)*0.5 ) &            ! d(ww)/dz
                             - ((pres0(i, j, k) - pres0(i, j, km))*dzhi(k))     ! - dp/dz
            end do
        end do
    end do
end subroutine advecw_Flimiter


subroutine computefluxlimiter(varp2,varp,var,varm,varm2,cL,cH,fluxL,fluxH)
    use modglobal, only: ilimiter, icentral, iupwind
    implicit none
    real, intent(in) :: varp2, varp, var, varm, varm2, cL, cH
    real, intent(out) :: fluxL, fluxH
    real :: r, fluxLlow, fluxLhigh, fluxHlow, fluxHhigh
    real, external :: fluxlimiter

    call computeflux(iupwind,varp,var,varm,cL,cH,fluxLlow,fluxHlow)
    call computeflux(icentral,varp,var,varm,cL,cH,fluxLhigh,fluxHhigh)
    
    if (cH .ge. 0.) then
        r = (var - varm) / (varp - var +1.0e-12)
    else
        r = (varp2 - varp) / (varp - var +1.0e-12)
    endif
    fluxH = fluxHlow - fluxlimiter(r,ilimiter) * (fluxHlow-fluxHhigh)

    if (cL .ge. 0.) then
        r = (varm - varm2) / (var - varm +1.0e-12)
    else
        r = (varp - var) / (var - varm +1.0e-12)
    end if
    fluxL = fluxLlow - fluxlimiter(r,ilimiter)*(fluxLlow-fluxLhigh);

end subroutine computefluxlimiter


subroutine computeflux(ischeme,varp,var,varm,cL,cH,fluxL,fluxH)
    use modglobal, only: icentral, iupwind
    implicit none
    integer, intent(in) :: ischeme
    real, intent(in) :: varp, var, varm, cL, cH
    real, intent(out) :: fluxL, fluxH
    select case (ischeme)
        case (icentral)      !! central difference
            fluxL = 0.5*cL*(varm+var)
            fluxH = 0.5*cH*(var+varp)
        case (iupwind)       !! upwind scheme
            if (cL .ge. 0.) then
                fluxL = cL*varm
            else
                fluxL = cL*var
            end if
            if (cH .ge. 0.) then
                fluxH = cH*var
            else
                fluxH = cH*varp
            end if
        case default
            write(0, *) "ERROR: Invalid discretization scheme"
            stop 1
    end select
end subroutine computeflux


subroutine compFLimiter_nonuni(varp2,varp,var,varm,varm2,delhp2,delhp,delh,delhm,delfp,delf,delfm,cL,cH,fluxL,fluxH)
    use modglobal, only: ilimiter, icentral, iupwind
    implicit none
    real, intent(in) :: varp2, varp, var, varm, varm2, delhp2, delhp, delh, delhm, delfp, delf, delfm, cL, cH
    real, intent(out) :: fluxL, fluxH
    real :: r, fluxLlow, fluxLhigh, fluxHlow, fluxHhigh
    real, external :: fluxlimiter

    call compF_nonuni(iupwind,varp,var,varm,delfp,delf,delfm,cL,cH,fluxLlow,fluxHlow)
    call compF_nonuni(icentral,varp,var,varm,delfp,delf,delfm,cL,cH,fluxLhigh,fluxHhigh)
    
    if (cH .ge. 0.) then
        r = ((var - varm)*delhp) / (delh*(varp - var +1.0e-12))
    else
        r = ((varp2 - varp)*delhp) / (delhp2*(varp - var +1.0e-12))
    endif
    fluxH = fluxHlow - fluxlimiter(r,ilimiter) * (fluxHlow-fluxHhigh)

    if (cL .ge. 0.) then
        r = ((varm - varm2)*delh) / (delhm*(var - varm +1.0e-12))
    else
        r = ((varp - var)*delh) / (delhp*(var - varm +1.0e-12))
    end if
    fluxL = fluxLlow - fluxlimiter(r,ilimiter)*(fluxLlow-fluxLhigh);

end subroutine compFLimiter_nonuni


subroutine compF_nonuni(ischeme,varp,var,varm,delp,del,delm,cL,cH,fluxL,fluxH)
    use modglobal, only: icentral, iupwind
    implicit none
    integer, intent(in) :: ischeme
    real, intent(in) :: varp, var, varm, delp, del, delm, cL, cH
    real, intent(out) :: fluxL, fluxH
    select case (ischeme)
        case (icentral)      !! central difference
            fluxL = cL*(varm*del+var*delm)/(del+delm)
            fluxH = cH*(var*delp+varp*del)/(delp+del)
        case (iupwind)       !! upwind scheme
            if (cL .ge. 0.) then
                fluxL = cL*varm
            else
                fluxL = cL*var
            end if
            if (cH .ge. 0.) then
                fluxH = cH*var
            else
                fluxH = cH*varp
            end if
        case default
            write(0, *) "ERROR: Invalid discretization scheme"
            stop 1
    end select
end subroutine compF_nonuni


real function fluxlimiter(r,ilimiterscheme)
    use modglobal, only: isuperbee, iminmod
    implicit none
    real, intent(in) :: r
    integer, intent(in) :: ilimiterscheme
    select case (ilimiterscheme)
        case (isuperbee)
            fluxlimiter = max( 0.0, max( min(1.0, 2.0*r), min(2.0, r) ) )
        case (iminmod)
            fluxlimiter = max( 0.0, min(1.0, r) )
        case default
            write(0, *) "ERROR: Unknown flux limiter scheme"
            stop 1
    end select
end function fluxlimiter