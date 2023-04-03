!  This file is part of uDALES.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 2006-2021 the uDALES Team.
!
SUBROUTINE wfuno(hi,hj,hk,iout1,iout2,iot,iomomflux,iotflux,iocth,obcTfluxA,utang1,utang2,Tcell,Twall,z0,z0h,n,ind,wforient)
   !wfuno
   !calculating wall function for momentum and scalars following Cai2012&Uno1995, extension of Louis 1979 method to rough walls
   !fluxes in m2/s2 and Km/s
   USE modglobal, ONLY : dzf,dzfi,dzh2i,dzhi,dzhiq,dy,dyi,dy2i,dyi5,dxf,dxh,dxfi,dxhi,dxh2i,ib,ie,jb,je,kb,ke,fkar,grav,jmax,rk3step,kmax,jge,jgb,dxfi5,dzfi5
   USE modsubgriddata, ONLY:ekh, ekm
   USE modmpi, ONLY:myid
   USE initfac, ONLY:block
   USE modibmdata
   use modfields, only : IIc, u0, v0, w0
   REAL, EXTERNAL :: unom
   INTEGER i, j, k, jl, ju, kl, ku, il, iu, km, im, jm, ip, jp, kp
   REAL :: Ribl0 = 0. !initial guess of Ribl based on Ts

   REAL :: bcTflux = 0. !temp storage for temperature flux
   REAL :: bcmomflux = 0. !temp storage for momentum flux
   REAL :: ctm = 0. !momentum transfer coefficient
   REAL :: cth = 0. !heat transfer coefficient
   REAL :: dummy = 0. !for debugging
   REAL :: delta = 0. !distance from wall
   REAL :: logdz = 0. !log(delta/z0)
   REAL :: logdzh = 0. !log(delta/z0h)
   REAL :: logzh = 0. !log(z0/z0h)
   REAL :: sqdz = 0. !sqrt(delta/z0)
   REAL :: utang1Int !Interpolated 1st tangential velocity component needed for stability calculation (to T location)
   REAL :: utang2Int !Interpolated 2nd tangential velocity component needed for stability calculation (to T location)
   REAL :: utangInt !Interpolated absolute tangential velocity
   REAL :: dT !Temperature difference between wall and cell
   REAL :: fkar2 !fkar^2, von Karman constant squared
   REAL :: emmo = 0., epmo = 0., epom = 0., emom = 0., eopm = 0., eomm = 0., empo = 0.
   REAL :: umin = 0.0001 !m^2/s^2

   INTEGER, INTENT(in) :: hi !<size of halo in i
   INTEGER, INTENT(in) :: hj !<size of halo in j
   INTEGER, INTENT(in) :: hk !<size of halo in k
   REAL, INTENT(out)   :: obcTfluxA !temperature flux of entire wall facet (double sum over indeces) [Km/s]
   REAL, INTENT(inout) :: iout1(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk) !updated prognostic tangential velocity (component1)
   REAL, INTENT(inout) :: iout2(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk) !updated prognostic tangential velocity (component2)
   REAL, INTENT(inout) :: iot(ib - hi:ie + hi, jb - hj:je + hj, kb:ke + hk) !updated prognostic temperature
   REAL, INTENT(inout) :: iomomflux(ib - hi:ie + hi, jb - hj:je + hj, kb-hk:ke + hk) !a field to save the momentum flux
   REAL, INTENT(inout) :: iotflux(ib - hi:ie + hi, jb - hj:je + hj, kb-hk:ke + hk) !a field to save the heat flux
   REAL, INTENT(inout) :: iocth(ib - hi:ie + hi, jb - hj:je + hj, kb-hk:ke + hk) !heat transfer coefficient, used to calculate moisture flux
   REAL, INTENT(in)    :: Tcell(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk) !Temperature of fluid cell
   REAL, INTENT(in)    :: Twall !Temperature of surfaces !SINCE EVERY WALL HAS PRECISELY ONE TEMPERATURE (at the outside). CAREFUL IF THIS EVER CHANGES (i.e. multiple EB facets per wall)
   REAL, INTENT(in)    :: z0
   REAL, INTENT(in)    :: z0h
   REAL, INTENT(in)    :: utang1(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk) !tangential velocity field
   REAL, INTENT(in)    :: utang2(ib - hi:ie + hi, jb - hj:je + hj, kb - hk:ke + hk) !second tangential velocity field
   INTEGER, INTENT(in) :: n ! number of the block, used to get i,j,k-indeces
   INTEGER, INTENT(in) :: ind ! in case of y-wall (case 3x & 4x) "ind" is used for j-index, otherwise this is irrelevant
   INTEGER, INTENT(in) :: wforient !orientation of the facet see below:
   !frist digit, orientation of wall, determines iteration indices
   !second digit, if for momentum or for scalar (necessary because of staggered grid -> which variable to interpolate)
   !xlow=1,xup=2,yup=3,ylow=4,z=5
   !momentum=1,scalar=2
   fkar2 = fkar**2
   obcTfluxA = 0.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CASES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CASES FOR SCALARS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SELECT CASE (wforient)
   CASE (12) !wall in yz -> wf in x (=vertical), lower wall, west wall
      i = block(n, 1) - 1 !wall property and fluid index
      ip = i + 1 !index to remove subgrid flux
      jl = MAX(block(n, 3) - myid*jmax, 1) ! starting j-index
      ju = MIN(block(n, 4) - myid*jmax, jmax) ! ending j-index
      kl = block(n, 5) ! starting k-index
      ku = block(n, 6) ! ending k-index
      delta = dxf(i)*0.5 !
      logdz = LOG(delta/z0)
      logdzh = LOG(delta/z0h)
      logzh = LOG(z0/z0h)
      sqdz = SQRT(delta/z0)

      DO k = kl, ku
         DO j = jl, ju
            utang1Int = (utang1(i, j, k) + utang1(i, j + 1, k))*0.5
            utang2Int = (utang2(i, j, k) + utang2(i, j, k + 1))*0.5
            utangInt = max((utang1Int**2 + utang2Int**2), umin)
            dT = (Tcell(i, j, k) - Twall)
            Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri

            call unoh(bcTflux, cth, logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
            !save result and update field
            obcTfluxA = obcTfluxA + bcTflux
            iocth(i,j,k) = cth
            iotflux(i, j, k) = iotflux(i, j, k) + bcTflux*dxfi(i)
            iot(i,j,k) = iot(i,j,k) - bcTflux*dxfi(i) &
                       - 0.5*(ekh(ip,j,k)*dxf(i) + ekh(i,j,k)*dxf(ip))*(Tcell(ip,j,k) - Tcell(i,j,k))*dxh2i(ip)*dxfi(i) & !
                       + (u0(ip, j, k)*(Tcell(ip, j, k)*dxf(i) + Tcell(i, j, k)*dxf(ip))*dxhi(ip))*dxfi5(i) &
                       - (u0(ip, j, k)*(Tcell(i , j, k)*dxf(i) + Tcell(i, j, k)*dxf(ip))*dxhi(ip))*dxfi5(i)
         END DO
      END DO

   CASE (22) ! wall in yz -> wf in x (=vertical), upper wall, east wall
      i = block(n, 2) + 1 !
      im = i - 1 !
      jl = MAX(block(n, 3) - myid*jmax, 1) ! starting j-index
      ju = MIN(block(n, 4) - myid*jmax, jmax) ! ending j-index
      kl = block(n, 5) ! starting k-index
      ku = block(n, 6) ! ending k-index
      delta = dxh(i)*0.5
      logdz = LOG(delta/z0)
      logdzh = LOG(delta/z0h)
      logzh = LOG(z0/z0h)
      sqdz = SQRT(delta/z0)

      DO k = kl, ku
         DO j = jl, ju
            utang1Int = (utang1(i, j, k) + utang1(i, j + 1, k))*0.5
            utang2Int = (utang2(i, j, k) + utang2(i, j, k + 1))*0.5
            utangInt = max(umin, (utang1Int**2 + utang2Int**2))
            dT = (Tcell(i, j, k) - Twall)
            Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri

            call unoh(bcTflux, cth, logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
            !save result and update field
            obcTfluxA = obcTfluxA + bcTflux
            iotflux(i, j, k) = iotflux(i, j, k) + bcTflux*dxfi(i)
            iocth(i,j,k) = cth
            iot(i,j,k) = iot(i,j,k) - bcTflux*dxfi(i) &
                       + 0.5*(ekh(i,j,k)*dxf(im) + ekh(im,j,k)*dxf(i))*(Tcell(i,j,k) - Tcell(im,j,k))*dxh2i(i)*dxfi(i) &
                       + (- u0(i, j, k)*(Tcell(im, j, k)*dxf(i) + Tcell(i, j, k)*dxf(im))*dxhi(i))*dxfi5(i) &
                       - (- u0(i, j, k)*(Tcell(i , j, k)*dxf(i) + Tcell(i, j, k)*dxf(im))*dxhi(i))*dxfi5(i)
         END DO
      END DO

   CASE (32) !wall in xz -> wf in y (=vertical) upper, north wall
      j = ind
      jm = j - 1
      il = block(n, 1)
      iu = block(n, 2)
      kl = block(n, 5)
      ku = block(n, 6)
      delta = 0.5*dy
      logdz = LOG(delta/z0)
      logdzh = LOG(delta/z0h)
      logzh = LOG(z0/z0h)
      sqdz = SQRT(delta/z0)

      DO k = kl, ku
         DO i = il, iu
            utang1Int = (utang1(i, j, k) + utang1(i + 1, j, k))*0.5
            utang2Int = (utang2(i, j, k) + utang2(i, j, k + 1))*0.5
            utangInt = max(umin, (utang1Int**2 + utang2Int**2))
            dT = (Tcell(i, j, k) - Twall)

            Ribl0 = grav*delta*dT/(Twall*utangInt) !

            call unoh(bcTflux, cth, logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
            obcTfluxA = obcTfluxA + bcTflux
            iotflux(i, j, k) = iotflux(i, j, k) + bcTflux*dyi
            iocth(i,j,k) = cth
            iot(i, j, k) = iot(i, j, k) - bcTflux*dyi &
                         + 0.5*(ekh(i, j, k) + ekh(i, jm, k))*(Tcell(i, j, k) - Tcell(i, jm, k))*dy2i &
                         + (- v0(i, j, k)*(Tcell(i, jm, k) + Tcell(i, j, k)))*dyi5 &
                         - (- v0(i, j, k)*(Tcell(i, j , k) + Tcell(i, j, k)))*dyi5
         END DO
      END DO

   CASE (42) !wall in xz -> wf in y (=vertical) lower, south wall
      j = ind
      jp = j + 1
      il = block(n, 1)
      iu = block(n, 2)
      kl = block(n, 5)
      ku = block(n, 6)
      delta = 0.5*dy
            logdz = LOG(delta/z0)
            logdzh = LOG(delta/z0h)
            logzh = LOG(z0/z0h)
            sqdz = SQRT(delta/z0)

      DO k = kl, ku
         DO i = il, iu
            utang1Int = (utang1(i, j, k) + utang1(i + 1, j, k))*0.5
            utang2Int = (utang2(i, j, k) + utang2(i, j, k + 1))*0.5
            utangInt = max(umin, (utang1Int**2 + utang2Int**2))
            dT = (Tcell(i, j, k) - Twall)

            Ribl0 = grav*delta*dT/(Twall*utangInt) !
            call unoh(bcTflux, cth, logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
            obcTfluxA = obcTfluxA + bcTflux
            iotflux(i, j, k) = iotflux(i, j, k) + bcTflux*dyi
            iocth(i,j,k) = cth
            iot(i, j, k) = iot(i, j, k) - bcTflux*dyi &
                         - 0.5*(ekh(i, jp, k) + ekh(i, j, k))*(Tcell(i, jp, k) - Tcell(i, j, k))*dy2i &
                         + v0(i, jp, k)*(Tcell(i, jp, k) + Tcell(i, j, k))*dyi5 &
                         - v0(i, jp, k)*(Tcell(i, j , k) + Tcell(i, j, k))*dyi5
         END DO
      END DO

   CASE (52) !wall in xy -> wf in z (=horizontal), top wall
      k = block(n, 6) + 1 !block location
      !if (.not.(k.gt.kmax)) then
      if (k > ke) return
      km = k - 1 !
      il = block(n, 1)
      iu = block(n, 2)
      jl = MAX(block(n, 3) - myid*jmax, 1)
      ju = MIN(block(n, 4) - myid*jmax, jmax)
      delta = dzf(k)*0.5
      logdz = LOG(delta/z0)
      logdzh = LOG(delta/z0h)
      logzh = LOG(z0/z0h)
      sqdz = SQRT(delta/z0)

      DO j = jl, ju
         DO i = il, iu
            utang1Int = (utang1(i, j, k) + utang1(i + 1, j, k))*0.5
            utang2Int = (utang2(i, j, k) + utang2(i, j + 1, k))*0.5
            utangInt = max(umin, (utang1Int**2 + utang2Int**2))
            dT = (Tcell(i, j, k) - Twall)

            Ribl0 = grav*delta*dT/(Twall*utangInt) !
            call unoh(bcTflux, cth, logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
            obcTfluxA = obcTfluxA + bcTflux
            iotflux(i, j, k) = iotflux(i, j, k) + bcTflux*dzfi(k)
            iocth(i,j,k) = cth
            iot(i, j, k) = iot(i, j, k) - bcTflux*dzfi(k) &
                         + 0.5*(dzf(km)*ekh(i, j, k) + dzf(k)*ekh(i, j, km))*(Tcell(i, j, k) - Tcell(i, j, km))*dzh2i(k)*dzfi(k) &
                         + (- w0(i, j, k)*(Tcell(i, j, km)*dzf(k) + Tcell(i, j, k)*dzf(km))*dzhi(k))*dzfi5(k) &
                         - (- w0(i, j, k)*(Tcell(i, j, k )*dzf(k) + Tcell(i, j, k)*dzf(km))*dzhi(k))*dzfi5(k)
         END DO
      END DO
      !end if

   ! CASE (62) !wall in xy -> wf in z (=horizontal), bottom wall
   !    k = block(n, 5) - 1 !block location - Should this have -1?
   !    !if (.not.(k.lt.0)) then
   !    if (k < kb) return
   !    kp = k + 1 !
   !    il = block(n, 1)
   !    iu = block(n, 2)
   !    jl = MAX(block(n, 3) - myid*jmax, 1)
   !    ju = MIN(block(n, 4) - myid*jmax, jmax)
   !    delta = dzf(k)*0.5
   !    logdz = LOG(delta/z0)
   !    logdzh = LOG(delta/z0h)
   !    logzh = LOG(z0/z0h)
   !    sqdz = SQRT(delta/z0)
   !
   !    DO j = jl, ju
   !      DO i = il, iu
   !        utang1Int = (utang1(i, j, k) + utang1(i + 1, j, k))*0.5
   !        utang2Int = (utang2(i, j, k) + utang2(i, j + 1, k))*0.5
   !        utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !        dT = (Tcell(i, j, k) - Twall)
   !
   !        Ribl0 = grav*delta*dT/(Twall*utangInt) !
   !        !call unoh(bcTflux, cth, logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
   !        obcTfluxA = obcTfluxA + bcTflux
   !        iotflux(i, j, k) = iotflux(i, j, k) + bcTflux*dzfi(k)
   !        iocth(i,j,k) = cth
   !        iot(i, j, k) = iot(i, j, k) - bcTflux*dzfi(k) &
   !                     - 0.5*(dzf(kp)*ekh(i, j, k) + dzf(k)*ekh(i, j, kp))*(Tcell(i, j, kp) - Tcell(i, j, k))*dzh2i(k)*dzfi(k) &
   !                     + w0(i, j, kp)*(Tcell(i, j, kp)*dzf(k) + Tcell(i, j, k)*dzf(kp))*dzhi(kp)*dzfi5(k) &
   !                     - w0(i, j, kp)*(Tcell(i, j, k )*dzf(k) + Tcell(i, j, k)*dzf(kp))*dzhi(kp)*dzfi5(k)
   !      END DO
   !    END DO
   !    !end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!CASES FOR MOMENTUM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   CASE (11) ! west wall
      i = block(n, 1) - 1 !fluid location (also where wall variables are stored)
      ip = i + 1 !inside wall, used for subtracting original diffusion term
      jl = MAX(block(n, 3) - myid*jmax, 1) + 1 ! starting j-index      !might cause problem when jl=1
      ju = MIN(block(n, 4) - myid*jmax, jmax) ! ending j-index     !might cause problem when ju=jmax
      kl = block(n, 5) ! starting k-index
      ku = block(n, 6) ! ending k-index

      delta = dxf(i)*0.5
      logdz = LOG(delta/z0)
      logdzh = LOG(delta/z0h)
      logzh = LOG(z0/z0h)
      sqdz = SQRT(delta/z0)

      !v west
      DO k = kl, ku
         DO j = jl, ju
            utang1Int = utang1(i, j, k)
            utang2Int = (utang2(i, j, k) + utang2(i, j, k + 1) + utang2(i, j - 1, k) + utang2(i, j - 1, k + 1))*0.25
            utangInt = max(umin, (utang1Int**2 + utang2Int**2))
            !dT = ((Tcell(i, j, k)*(2-IIc(i, j - 1, k) + Tcell(i, j - 1, k)*) - (Twall + Twall))*0.5
            dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j - 1, k)) + Tcell(i, j - 1, k)*IIc(i, j - 1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
            Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
            ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
            !dummy = (utang1Int**2)*ctm
            dummy = abs(utang1Int)*sqrt(utangInt)*ctm
            bcmomflux = SIGN(dummy, utang1Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)

            epmo = 0.25*((ekm(i, j, k) + ekm(i, j - 1, k))*dxf(ip) + &
                         (ekm(ip, j, k) + ekm(ip, j - 1, k))*dxf(i))*dxhi(ip)

            iout1(i, j, k) = iout1(i, j, k) - (utang1(ip, j, k) - utang1(i, j, k))*epmo*dxhi(ip)*dxfi(i) - bcmomflux*dxfi(i) !

         END DO
      END DO

      !v west edge south
      j = MAX(block(n, 3) - myid*jmax, 1)
      DO k = kl, ku
         utang1Int = utang1(i, j, k)
         utang2Int = (utang2(i, j, k) + utang2(i, j, k + 1) + utang2(i, j - 1, k) + utang2(i, j - 1, k + 1))*0.25
         utangInt = max(umin, (utang1Int**2 + utang2Int**2))
         !dT = Tcell(i, j, k) - Twall
         dT = (Tcell(i, j, k) - Twall)*IIc(i,j,k)
         Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri
         !dummy = (utang1Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
         ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
         dummy = abs(utang1Int)*sqrt(utangInt)*ctm
         bcmomflux = SIGN(dummy, utang1Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)*0.5 !only half since on edge of block (another half might come from another processor?)

         !   epmo = 0.5 * (ekm(ip,j,k)*dxf(i) + ekm(i,j,k)*dxf(ip)) * dxhi(ip)
         epmo = 0.25*((ekm(i, j, k) + ekm(i, j - 1, k))*dxf(ip) + &
                      (ekm(ip, j, k) + ekm(ip, j - 1, k))*dxf(i))*dxhi(ip)
         iout1(i, j, k) = iout1(i, j, k) - ((utang1(ip, j, k) - utang1(i, j, k))*epmo*dxhi(ip)*dxfi(i) - bcmomflux*dxfi(i))*0.5 ! remove standard diffusion apply only half of wall-flux since it's an edge
         !only half of the flux, since only half of the control-volume around v is touching this facet (other half is either in free air or touching another facet)
      END DO

      !v west edge north
      j = MIN(block(n, 4) - myid*jmax, jmax) + 1
      DO k = kl, ku
         utang1Int = utang1(i, j, k)
         utang2Int = (utang2(i, j, k) + utang2(i, j, k + 1) + utang2(i, j - 1, k) + utang2(i, j - 1, k + 1))*0.25
         utangInt = max(umin, (utang1Int**2 + utang2Int**2))
         !dT = Tcell(i, j - 1, k) - Twall !possibly look at this
         dT = (Tcell(i,j-1,k)-Twall)*IIc(i,j-1,k)
         Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri
         !dummy = (utang1Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
         ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
         dummy = abs(utang1Int)*sqrt(utangInt)*ctm
         bcmomflux = SIGN(dummy, utang1Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)*0.5
         epmo = 0.25*((ekm(i, j, k) + ekm(i, j - 1, k))*dxf(ip) + &
                      (ekm(ip, j, k) + ekm(ip, j - 1, k))*dxf(i))*dxhi(ip)
         iout1(i, j, k) = iout1(i, j, k) - (utang1(ip, j, k) - utang1(i, j, k))*epmo*dxhi(ip)*dxfi(i) - bcmomflux*dxfi(i)*0.5 ! %remove standard diffusion apply only half of wall-flux since it's an edge

      END DO

      !w west
      jl = MAX(block(n, 3) - myid*jmax, 1) !
      kl = block(n, 5) + 1 !
      DO k = kl, ku
         DO j = jl, ju
            utang1Int = (utang1(i, j, k) + utang1(i, j + 1, k) + utang1(i, j + 1, k - 1) + utang1(i, j, k - 1))*0.25
            utang2Int = utang2(i, j, k)
            utangInt = max(umin, (utang1Int**2 + utang2Int**2))
            !dT = ((Tcell(i, j, k) + Tcell(i, j, k - 1)) - (Twall + Twall))*0.5
            dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j, k-1)) + Tcell(i, j, k-1)*IIc(i, j, k-1)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
            Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
            ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
            !dummy = (utang2Int**2)*ctm
            dummy = abs(utang2Int)*sqrt(utangInt)*ctm
            bcmomflux = SIGN(dummy, utang2Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)
            epom = (dzf(k - 1)*(ekm(i, j, k)*dxf(ip) + ekm(ip, j, k)*dxf(i))*dxhi(ip) + &
                    dzf(k)*(ekm(i, j, k - 1)*dxf(ip) + ekm(ip, j, k - 1)*dxf(i))*dxhi(ip))*dzhiq(k)
            iout2(i, j, k) = iout2(i, j, k) - (utang2(ip, j, k) - utang2(i, j, k))*epom*dxhi(ip)*dxfi(i) - bcmomflux*dxfi(i) !
         END DO
      END DO

      !w west top edge
      k = block(n, 6) + 1 ! ending k-index
      km = k - 1
      DO j = jl, ju
         utang1Int = (utang1(i, j, k) + utang1(i, j + 1, k) + utang1(i, j + 1, k - 1) + utang1(i, j, k - 1))*0.25
         utang2Int = utang2(i, j, k)
         utangInt = max(umin, (utang1Int**2 + utang2Int**2))
         !dT = Tcell(i, j, k - 1) - Twall
         dT = (Tcell(i,j,k-1)-Twall)*IIc(i,j,k-1)

         Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri
         !dummy = (utang2Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
         ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
         dummy = abs(utang2Int)*sqrt(utangInt)*ctm
         bcmomflux = SIGN(dummy, utang2Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)*0.5
         epom = (dzf(km)*(ekm(i, j, k)*dxf(ip) + ekm(ip, j, k)*dxf(i))*dxhi(ip) + &
                 dzf(k)*(ekm(i, j, km)*dxf(ip) + ekm(ip, j, km)*dxf(i))*dxhi(ip))*dzhiq(k)
         iout2(i, j, k) = iout2(i, j, k) - (utang2(ip, j, k) - utang2(i, j, k))*epom*dxhi(ip)*dxfi(i) - bcmomflux*dxfi(i)*0.5 !
      END DO

      !w west bottom edge
      k = block(n, 5)  ! ending k-index
      !if (k.gt.0) then
      if (k > kb) then
      km = k - 1
      DO j = jl, ju
         utang1Int = (utang1(i, j, k) + utang1(i, j + 1, k) + utang1(i, j + 1, k - 1) + utang1(i, j, k - 1))*0.25
         utang2Int = utang2(i, j, k)
         utangInt = max(umin, (utang1Int**2 + utang2Int**2))
         !dT = Tcell(i, j, k - 1) - Twall
         dT = (Tcell(i,j,k-1)-Twall)*IIc(i,j,k-1)
         Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri
         !dummy = (utang2Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
         ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
         dummy = abs(utang2Int)*sqrt(utangInt)*ctm
         bcmomflux = SIGN(dummy, utang2Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)*0.5
         epom = (dzf(km)*(ekm(i, j, k)*dxf(ip) + ekm(ip, j, k)*dxf(i))*dxhi(ip) + &
                 dzf(k)*(ekm(i, j, km)*dxf(ip) + ekm(ip, j, km)*dxf(i))*dxhi(ip))*dzhiq(k)
         iout2(i, j, k) = iout2(i, j, k) - (utang2(ip, j, k) - utang2(i, j, k))*epom*dxhi(ip)*dxfi(i) - bcmomflux*dxfi(i)*0.5 !
      END DO
      end if

   CASE (21)
      !v east
      i = block(n, 2) + 1 !fluid
      im = i - 1 !inside block
      jl = MAX(block(n, 3) - myid*jmax, 1) + 1 ! starting j-index      !might cause problem when jl=1
      ju = MIN(block(n, 4) - myid*jmax, jmax) ! ending j-index     !might cause problem when ju=jmax
      kl = block(n, 5) ! starting k-index
      ku = block(n, 6) ! ending k-index
      delta = dxh(i)*0.5
      logdz = LOG(delta/z0)
      logdzh = LOG(delta/z0h)
      logzh = LOG(z0/z0h)
      sqdz = SQRT(delta/z0)

      DO k = kl, ku
         DO j = jl, ju
            utang1Int = utang1(i, j, k)
            utang2Int = (utang2(i, j, k) + utang2(i, j, k + 1) + utang2(i, j - 1, k) + utang2(i, j - 1, k + 1))*0.25
            utangInt = max(umin, (utang1Int**2 + utang2Int**2))
            !dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
            dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j - 1, k)) + Tcell(i, j - 1, k)*IIc(i, j - 1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
            Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
            !call function repeatedly
            ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
            !dummy = (utang1Int**2)*ctm
            dummy = abs(utang1Int)*sqrt(utangInt)*ctm
            bcmomflux = SIGN(dummy, utang1Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)
            emmo = 0.25*((ekm(i, j, k) + ekm(i, j - 1, k))*dxf(im) + (ekm(im, j - 1, k) + ekm(im, j, k))*dxf(i))*dxhi(i) ! dx is non-equidistant
            iout1(i, j, k) = iout1(i, j, k) + (utang1(i, j, k) - utang1(im, j, k))*emmo*dxhi(i)*dxfi(i) - bcmomflux*dxfi(i) !
         END DO
      END DO

      !v east edge south
      j = MAX(block(n, 3) - myid*jmax, 1)
      DO k = kl, ku
         utang1Int = utang1(i, j, k)
         utang2Int = (utang2(i, j, k) + utang2(i, j, k + 1) + utang2(i, j - 1, k) + utang2(i, j - 1, k + 1))*0.25
         utangInt = max(umin, (utang1Int**2 + utang2Int**2))
         !dT = Tcell(i, j, k) - Twall
         dT = (Tcell(i, j, k) - Twall)*IIc(i,j,k)
         Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri
         !dummy = (utang1Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
         ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
         dummy = abs(utang1Int)*sqrt(utangInt)*ctm
         bcmomflux = SIGN(dummy, utang1Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)*0.5
         emmo = 0.25*((ekm(i, j, k) + ekm(i, j - 1, k))*dxf(im) + (ekm(im, j - 1, k) + ekm(im, j, k))*dxf(i))*dxhi(i) ! dx is non-equidistant
         iout1(i, j, k) = iout1(i, j, k) + (utang1(i, j, k) - utang1(im, j, k))*emmo*dxhi(i)*dxfi(i) - bcmomflux*dxfi(i)*0.5 !
      END DO

      !v east edge north
      j = MIN(block(n, 4) - myid*jmax, jmax) + 1 !
      DO k = kl, ku
         utang1Int = utang1(i, j, k)
         utang2Int = (utang2(i, j, k) + utang2(i, j, k + 1) + utang2(i, j - 1, k) + utang2(i, j - 1, k + 1))*0.25
         utangInt = max(umin, (utang1Int**2 + utang2Int**2))
         !dT = Tcell(i, j - 1, k) - Twall
         dT = (Tcell(i, j-1, k) - Twall)*IIc(i, j-1, k)
         Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri
         !dummy = (utang1Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
         ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
         dummy = abs(utang1Int)*sqrt(utangInt)*ctm
         bcmomflux = SIGN(dummy, utang1Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)*0.5
         emmo = 0.25*((ekm(i, j, k) + ekm(i, j - 1, k))*dxf(im) + (ekm(im, j - 1, k) + ekm(im, j, k))*dxf(i))*dxhi(i) ! dx is non-equidistant
         iout1(i, j, k) = iout1(i, j, k) + (utang1(i, j, k) - utang1(im, j, k))*emmo*dxhi(i)*dxfi(i) - bcmomflux*dxfi(i)*0.5 !
      END DO

      !w east
      jl = MAX(block(n, 3) - myid*jmax, 1) !
      kl = block(n, 5) + 1 !
      DO k = kl, ku
         DO j = jl, ju
            utang1Int = (utang1(i, j, k) + utang1(i, j + 1, k) + utang1(i, j + 1, k - 1) + utang1(i, j, k - 1))*0.25
            utang2Int = utang2(i, j, k)
            utangInt = max(umin, (utang1Int**2 + utang2Int**2))
            !dT = ((Tcell(i, j, k) + Tcell(i, j, k - 1)) - (Twall + Twall))*0.5
            dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j, k-1)) + Tcell(i, j, k-1)*IIc(i, j, k-1)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
            Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
            ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
            !dummy = (utang2Int**2)*ctm
            dummy = abs(utang2Int)*sqrt(utangInt)*ctm
            bcmomflux = SIGN(dummy, utang2Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)
            emom = (dzf(k - 1)*(ekm(i, j, k)*dxf(im) + ekm(im, j, k)*dxf(i))*dxhi(i) + &
                    dzf(k)*(ekm(i, j, k - 1)*dxf(im) + ekm(im, j, k - 1)*dxf(i))*dxhi(i))*dzhiq(k)
            iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(im, j, k))*emom*dxhi(i)*dxfi(i) - bcmomflux*dxfi(i) !
         END DO
      END DO

      !w east edge top
      k = block(n, 6) + 1 ! ending k-index
      DO j = jl, ju
         utang1Int = (utang1(i, j, k) + utang1(i, j + 1, k) + utang1(i, j + 1, k - 1) + utang1(i, j, k - 1))*0.25
         utang2Int = utang2(i, j, k)
         utangInt = max(umin, (utang1Int**2 + utang2Int**2))
         !dT = Tcell(i, j, k - 1) - Twall
         dT = (Tcell(i,j,k-1)-Twall)*IIc(i,j,k-1)
         Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri
         !dummy = (utang2Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
         ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
         dummy = abs(utang2Int)*sqrt(utangInt)*ctm
         bcmomflux = SIGN(dummy, utang2Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)*0.5
         emom = (dzf(k - 1)*(ekm(i, j, k)*dxf(im) + ekm(im, j, k)*dxf(i))*dxhi(i) + &
                 dzf(k)*(ekm(i, j, k - 1)*dxf(im) + ekm(im, j, k - 1)*dxf(i))*dxhi(i))*dzhiq(k)
         iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(im, j, k))*emom*dxhi(i)*dxfi(i) - bcmomflux*dxfi(i)*0.5 !
      END DO

     !w east edge bot
      k = block(n, 5)  !
      !if (k.gt.0) then
      if (k > kb) then
      DO j = jl, ju
         utang1Int = (utang1(i, j, k) + utang1(i, j + 1, k) + utang1(i, j + 1, k - 1) + utang1(i, j, k - 1))*0.25
         utang2Int = utang2(i, j, k)
         utangInt = max(umin, (utang1Int**2 + utang2Int**2))
         !dT = Tcell(i, j, k - 1) - Twall
         dT = (Tcell(i,j,k-1)-Twall)*IIc(i,j,k-1)
         Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri
         !dummy = (utang2Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
         ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
         dummy = abs(utang2Int)*sqrt(utangInt)*ctm
         bcmomflux = SIGN(dummy, utang2Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dxfi(i)*0.5
         emom = (dzf(k - 1)*(ekm(i, j, k)*dxf(im) + ekm(im, j, k)*dxf(i))*dxhi(i) + &
         dzf(k)*(ekm(i, j, k - 1)*dxf(im) + ekm(im, j, k - 1)*dxf(i))*dxhi(i))*dzhiq(k)
         iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(im, j, k))*emom*dxhi(i)*dxfi(i) - bcmomflux*dxfi(i)*0.5 !
      END DO
      end if

   CASE (31) !wall in xz -> wf in y (=vertical) upper, north wall
      j = ind
      jm = j - 1
      il = block(n, 1) + 1
      iu = block(n, 2)
      kl = block(n, 5)
      ku = block(n, 6)
      delta = 0.5*dy
      logdz = LOG(delta/z0)
      logdzh = LOG(delta/z0h)
      logzh = LOG(z0/z0h)
      sqdz = SQRT(delta/z0)

      !u north
      DO k = kl, ku
         DO i = il, iu
            utang1Int = utang1(i, j, k)
            utang2Int = (utang2(i, j, k) + utang2(i - 1, j, k) + utang2(i, j, k + 1) + utang2(i - 1, j, k + 1))*0.25
            utangInt = max(umin, (utang1Int**2 + utang2Int**2))
            !dT = ((Tcell(i, j, k) + Tcell(i - 1, j, k)) - (Twall + Twall))*0.5
            dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i-1, j, k)) + Tcell(i-1, j, k)*IIc(i-1, j, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
            Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
            ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
            !dummy = (utang1Int**2)*ctm
            dummy = abs(utang1Int)*sqrt(utangInt)*ctm
            bcmomflux = SIGN(dummy, utang1Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi

            emmo = 0.25*((ekm(i, j, k) + ekm(i, jm, k))*dxf(i - 1) + (ekm(i - 1, jm, k) + ekm(i - 1, j, k))*dxf(i))*dxhi(i)

            iout1(i, j, k) = iout1(i, j, k) + (utang1(i, j, k) - utang1(i, jm, k))*emmo*dy2i-bcmomflux*dyi !
         END DO
      END DO

      !u north east edge
      i = block(n, 2) + 1
      DO k = kl, ku
         utang1Int = utang1(i, j, k)
         utang2Int = (utang2(i, j, k) + utang2(i - 1, j, k) + utang2(i, j, k + 1) + utang2(i - 1, j, k + 1))*0.25
         utangInt = max(umin, (utang1Int**2 + utang2Int**2))
         !dT = Tcell(i - 1, j, k) - Twall
         dT = (Tcell(i-1,j,k)-Twall)*IIc(i-1,j,k)
         Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri
         !dummy = (utang1Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
         ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
         dummy = abs(utang1Int)*sqrt(utangInt)*ctm
         bcmomflux = SIGN(dummy, utang1Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi5
         emmo = 0.25*((ekm(i, j, k) + ekm(i, jm, k))*dxf(i - 1) + (ekm(i - 1, jm, k) + ekm(i - 1, j, k))*dxf(i))*dxhi(i)
         iout1(i, j, k) = iout1(i, j, k) + (utang1(i, j, k) - utang1(i, jm, k))*emmo*dy2i-bcmomflux*dyi*0.5 !
      END DO

      !u north west edge
      i = block(n, 1)
      DO k = kl, ku
         utang1Int = utang1(i, j, k)
         utang2Int = (utang2(i, j, k) + utang2(i - 1, j, k) + utang2(i, j, k + 1) + utang2(i - 1, j, k + 1))*0.25
         utangInt = max(umin, (utang1Int**2 + utang2Int**2))
         !dT = Tcell(i, j, k) - Twall
         dT = (Tcell(i,j,k)-Twall)*IIc(i,j,k)
         Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri
         !dummy = (utang1Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
         ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
         dummy = abs(utang1Int)*sqrt(utangInt)*ctm
         bcmomflux = SIGN(dummy, utang1Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi5
         emmo = 0.25*((ekm(i, j, k) + ekm(i, jm, k))*dxf(i - 1) + (ekm(i - 1, jm, k) + ekm(i - 1, j, k))*dxf(i))*dxhi(i)
         iout1(i, j, k) = iout1(i, j, k) + (utang1(i, j, k) - utang1(i, jm, k))*emmo*dy2i-bcmomflux*dyi*0.5 !
      END DO

      !w north
      il = block(n, 1) !
      kl = block(n, 5) + 1 !
      DO k = kl, ku
         DO i = il, iu
            utang1Int = (utang1(i, j, k) + utang1(i, j, k - 1) + utang1(i + 1, j, k) + utang1(i + 1, j, k - 1))*0.25
            utang2Int = utang2(i, j, k)
            utangInt = max(umin, (utang1Int**2 + utang2Int**2))
            !dT = ((Tcell(i, j, k) + Tcell(i, j, k - 1)) - (Twall + Twall))*0.5
            dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j, k-1)) + Tcell(i, j, k-1)*IIc(i, j, k-1)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
            Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
            ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
            !dummy = (utang2Int**2)*ctm
            dummy = abs(utang2Int)*sqrt(utangInt)*ctm
            bcmomflux = SIGN(dummy, utang2Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi
            eomm = (dzf(k - 1)*(ekm(i, j, k) + ekm(i, jm, k)) + dzf(k)*(ekm(i, j, k - 1) + ekm(i, jm, k - 1)))*dzhiq(k) ! dz is non-eqidistant
            iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, jm, k))*eomm*dy2i-bcmomflux*dyi !
         END DO
      END DO

      !w north edge top
      k = block(n, 6) + 1
      DO i = il, iu
         utang1Int = (utang1(i, j, k) + utang1(i, j, k - 1) + utang1(i + 1, j, k) + utang1(i + 1, j, k - 1))*0.25
         utang2Int = utang2(i, j, k)
         utangInt = max(umin, (utang1Int**2 + utang2Int**2))
         !dT = Tcell(i, j, k - 1) - Twall
         dT = (Tcell(i,j,k-1)-Twall)*IIc(i,j,k-1)
         Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri
         !dummy = (utang2Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
         ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
         dummy = abs(utang2Int)*sqrt(utangInt)*ctm
         bcmomflux = SIGN(dummy, utang2Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi5
         eomm = (dzf(k - 1)*(ekm(i, j, k) + ekm(i, jm, k)) + dzf(k)*(ekm(i, j, k - 1) + ekm(i, jm, k - 1)))*dzhiq(k) ! dz is non-eqidistant
         iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, jm, k))*eomm*dy2i-bcmomflux*dyi*0.5 !
      END DO

      !w north edge bot
      k = block(n, 5)
      !if (k.gt.0) then
      if (k > kb) then
      DO i = il, iu
         utang1Int = (utang1(i, j, k) + utang1(i, j, k - 1) + utang1(i + 1, j, k) + utang1(i + 1, j, k - 1))*0.25
         utang2Int = utang2(i, j, k)
         utangInt = max(umin, (utang1Int**2 + utang2Int**2))
         !dT = Tcell(i, j, k - 1) - Twall
         dT = (Tcell(i, j, k - 1) - Twall)*IIc(i,j,k-1)
         Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri
         !dummy = (utang2Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
         ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
         dummy = abs(utang2Int)*sqrt(utangInt)*ctm
         bcmomflux = SIGN(dummy, utang2Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi5
         eomm = (dzf(k - 1)*(ekm(i, j, k) + ekm(i, jm, k)) + dzf(k)*(ekm(i, j, k - 1) + ekm(i, jm, k - 1)))*dzhiq(k) ! dz is non-eqidistant
         iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, jm, k))*eomm*dy2i-bcmomflux*dyi*0.5 !
      END DO
      end if

   CASE (41) !wall in xz -> wf in y (=vertical) lower, south wall
      j = ind
      jp = j + 1
      il = block(n, 1) + 1
      iu = block(n, 2)
      kl = block(n, 5)
      ku = block(n, 6)
      delta = 0.5*dy
      logdz = LOG(delta/z0)
      logdzh = LOG(delta/z0h)
      logzh = LOG(z0/z0h)
      sqdz = SQRT(delta/z0)

      DO k = kl, ku
         DO i = il, iu
            utang1Int = utang1(i, j, k)
            utang2Int = (utang2(i, j, k) + utang2(i - 1, j, k) + utang2(i, j, k + 1) + utang2(i - 1, j, k + 1))*0.25
            utangInt = max(umin, (utang1Int**2 + utang2Int**2))
            !dT = ((Tcell(i, j, k) + Tcell(i - 1, j, k)) - (Twall + Twall))*0.5
            dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i-1, j, k)) + Tcell(i-1, j, k)*IIc(i-1, j, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
            Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
            ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
            !dummy = (utang1Int**2)*ctm
            dummy = abs(utang1Int)*sqrt(utangInt)*ctm
            bcmomflux = SIGN(dummy, utang1Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi
            empo = 0.25*((ekm(i, j, k) + ekm(i, jp, k))*dxf(i - 1) + &
                         (ekm(i - 1, j, k) + ekm(i - 1, jp, k))*dxf(i))*dxhi(i) ! dx is non-equidistant
            iout1(i, j, k) = iout1(i, j, k) - (utang1(i, jp, k) - utang1(i, j, k))*empo*dy2i-bcmomflux*dyi !
         END DO
      END DO

      !u south edge west
      i = block(n, 1)
      DO k = kl, ku
         utang1Int = utang1(i, j, k)
         utang2Int = (utang2(i, j, k) + utang2(i - 1, j, k) + utang2(i, j, k + 1) + utang2(i - 1, j, k + 1))*0.25
         utangInt = max(umin, (utang1Int**2 + utang2Int**2))
         !dT = Tcell(i, j, k) - Twall
         dT = (Tcell(i,j,k)-Twall)*IIc(i,j,k)
         Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri
         !dummy = (utang1Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
         ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
         dummy = abs(utang1Int)*sqrt(utangInt)*ctm
         bcmomflux = SIGN(dummy, utang1Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi5
         empo = 0.25*((ekm(i, j, k) + ekm(i, jp, k))*dxf(i - 1) + &
                      (ekm(i - 1, j, k) + ekm(i - 1, jp, k))*dxf(i))*dxhi(i) ! dx is non-equidistant
         iout1(i, j, k) = iout1(i, j, k) - (utang1(i, jp, k) - utang1(i, j, k))*empo*dy2i-bcmomflux*dyi*0.5 !
      END DO

      !u south edge east
      i = block(n, 2) + 1
      DO k = kl, ku
         utang1Int = utang1(i, j, k)
         utang2Int = (utang2(i, j, k) + utang2(i - 1, j, k) + utang2(i, j, k + 1) + utang2(i - 1, j, k + 1))*0.25
         utangInt = max(umin, (utang1Int**2 + utang2Int**2))
         !dT = Tcell(i - 1, j, k) - Twall
         dT = (Tcell(i-1,j,k)-Twall)*IIc(i-1,j,k)
         Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri
         !dummy = (utang1Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
         ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
         dummy = abs(utang1Int)*sqrt(utangInt)*ctm
         bcmomflux = SIGN(dummy, utang1Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi5
         empo = 0.25*((ekm(i, j, k) + ekm(i, jp, k))*dxf(i - 1) + &
                      (ekm(i - 1, j, k) + ekm(i - 1, jp, k))*dxf(i))*dxhi(i) ! dx is non-equidistant
         iout1(i, j, k) = iout1(i, j, k) - (utang1(i, jp, k) - utang1(i, j, k))*empo*dy2i-bcmomflux*dyi*0.5 !
      END DO

      !w south
      il = block(n, 1) !
      kl = block(n, 5) + 1 !
      DO k = kl, ku
         DO i = il, iu
            utang1Int = (utang1(i, j, k) + utang1(i, j, k - 1) + utang1(i + 1, j, k) + utang1(i + 1, j, k - 1))*0.25
            utang2Int = utang2(i, j, k)
            utangInt = max(umin, (utang1Int**2 + utang2Int**2))
            !dT = ((Tcell(i, j, k) + Tcell(i, j, k - 1)) - (Twall + Twall))*0.5
            dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j, k-1)) + Tcell(i, j, k-1)*IIc(i, j, k-1)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
            Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
            ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
            !dummy = (utang2Int**2)*ctm
            dummy = abs(utang2Int)*sqrt(utangInt)*ctm
            bcmomflux = SIGN(dummy, utang2Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi
            eopm = (dzf(k - 1)*(ekm(i, j, k) + ekm(i, jp, k)) + &
                    dzf(k)*(ekm(i, j, k - 1) + ekm(i, jp, k - 1)))*dzhiq(k)
            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, jp, k) - utang2(i, j, k))*eopm*dy2i-bcmomflux*dyi !
         END DO
      END DO

      !w south edge top
      k = block(n, 6) + 1
      DO i = il, iu
         utang1Int = (utang1(i, j, k) + utang1(i, j, k - 1) + utang1(i + 1, j, k) + utang1(i + 1, j, k - 1))*0.25
         utang2Int = utang2(i, j, k)
         utangInt = max(umin, (utang1Int**2 + utang2Int**2))
         !dT = Tcell(i, j, k - 1) - Twall
         dT = (Tcell(i,j,k-1)-Twall)*IIc(i,j,k-1)
         Ribl0 = grav*delta*dT*2/(Twall*utangInt) !Eq. 6, guess initial Ri
         !call function repeatedly
         !dummy = (utang2Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
         ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
         dummy = abs(utang2Int)*sqrt(utangInt)*ctm
         bcmomflux = SIGN(dummy, utang2Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi5
         eopm = (dzf(k - 1)*(ekm(i, j, k) + ekm(i, jp, k)) + &
                 dzf(k)*(ekm(i, j, k - 1) + ekm(i, jp, k - 1)))*dzhiq(k)

         iout2(i, j, k) = iout2(i, j, k) - (utang2(i, jp, k) - utang2(i, j, k))*eopm*dy2i-bcmomflux*dyi*0.5 !
      END DO

      !w south edge bot
      k = block(n, 5)
      !if (k.gt.0) then
      if (k > kb) then
      DO i = il, iu
         utang1Int = (utang1(i, j, k) + utang1(i, j, k - 1) + utang1(i + 1, j, k) + utang1(i + 1, j, k - 1))*0.25
         utang2Int = utang2(i, j, k)
         utangInt = max(umin, (utang1Int**2 + utang2Int**2))
         !dT = Tcell(i, j, k - 1) - Twall
         dT = (Tcell(i,j,k-1)-Twall)*IIc(i,j,k-1)
         Ribl0 = grav*delta*dT*2/(Twall*utangInt) !Eq. 6, guess initial Ri
         !dummy = (utang2Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
         ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
         dummy = abs(utang2Int)*sqrt(utangInt)*ctm
         bcmomflux = SIGN(dummy, utang2Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dyi5
         eopm = (dzf(k - 1)*(ekm(i, j, k) + ekm(i, jp, k)) + &
         dzf(k)*(ekm(i, j, k - 1) + ekm(i, jp, k - 1)))*dzhiq(k)
         iout2(i, j, k) = iout2(i, j, k) - (utang2(i, jp, k) - utang2(i, j, k))*eopm*dy2i-bcmomflux*dyi*0.5 !
      END DO
      end if

   CASE (51) !wall in xy -> wf in z (=horizontal), top wall
      k = block(n, 6) + 1 !block location
      !if (.not.(k.gt.kmax)) then
      if (k > ke) return
      km = k - 1 !shear velocity location
      il = block(n, 1) + 1
      iu = block(n, 2)
      jl = MAX(block(n, 3) - myid*jmax, 1)
      ju = MIN(block(n, 4) - myid*jmax, jmax)
      delta = 0.5*dzf(k)
      logdz = LOG(delta/z0)
      logdzh = LOG(delta/z0h)
      logzh = LOG(z0/z0h)
      sqdz = SQRT(delta/z0)

      DO j = jl, ju
         DO i = il, iu
            utang1Int = utang1(i, j, k)
            utang2Int = (utang2(i, j, k) + utang2(i - 1, j, k) + utang2(i, j + 1, k) + utang2(i - 1, j + 1, k))*0.25
            utangInt = max(umin, (utang1Int**2 + utang2Int**2))
            !dT = ((Tcell(i, j, k) + Tcell(i - 1, j, k)) - (Twall + Twall))*0.5
            dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i-1, j, k)) + Tcell(i-1, j, k)*IIc(i-1, j, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
            Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
            ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
            !dummy = (utang1Int**2)*ctm
            dummy = abs(utang1Int)*sqrt(utangInt)*ctm
            bcmomflux = SIGN(dummy, utang1Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
            emom = (dzf(km)*(ekm(i, j, k)*dxf(i - 1) + ekm(i - 1, j, k)*dxf(i)) + &
            dzf(k)*(ekm(i, j, km)*dxf(i - 1) + ekm(i - 1, j, km)*dxf(i)))*dxhi(i)*dzhiq(k)
            iout1(i, j, k) = iout1(i, j, k) + (utang1(i, j, k) - utang1(i, j, km))*emom*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k) !
         END DO
      END DO

      !u top edge west
      i = block(n, 1)
      DO j = jl, ju
         utang1Int = utang1(i, j, k)
         utang2Int = (utang2(i, j, k) + utang2(i - 1, j, k) + utang2(i, j + 1, k) + utang2(i - 1, j + 1, k))*0.25
         utangInt = max(umin, (utang1Int**2 + utang2Int**2))
         !dT = Tcell(i, j, k) - Twall
         dT = (Tcell(i,j,k)-Twall)*IIc(i,j,k)
         Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri
         !dummy = (utang1Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
         ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
         dummy = abs(utang1Int)*sqrt(utangInt)*ctm
         bcmomflux = SIGN(dummy, utang1Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)*0.5
         emom = (dzf(km)*(ekm(i, j, k)*dxf(i - 1) + ekm(i - 1, j, k)*dxf(i)) + &
         dzf(k)*(ekm(i, j, km)*dxf(i - 1) + ekm(i - 1, j, km)*dxf(i)))*dxhi(i)*dzhiq(k)
         iout1(i, j, k) = iout1(i, j, k) + (utang1(i, j, k) - utang1(i, j, km))*emom*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5 !
      END DO

       !u top edge east
       i = block(n, 2) + 1
       DO j = jl, ju
         utang1Int = utang1(i, j, k)
         utang2Int = (utang2(i, j, k) + utang2(i - 1, j, k) + utang2(i, j + 1, k) + utang2(i - 1, j + 1, k))*0.25
         utangInt = max(umin, (utang1Int**2 + utang2Int**2))
         !dT = Tcell(i - 1, j, k) - Twall
         dT = (Tcell(i-1,j,k)-Twall)*IIc(i-1,j,k)
         Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri
         !dummy = (utang1Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
         ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
         dummy = abs(utang1Int)*sqrt(utangInt)*ctm
         bcmomflux = SIGN(dummy, utang1Int)
         iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)*0.5
         emom = (dzf(km)*(ekm(i, j, k)*dxf(i - 1) + ekm(i - 1, j, k)*dxf(i)) + &
         dzf(k)*(ekm(i, j, km)*dxf(i - 1) + ekm(i - 1, j, km)*dxf(i)))*dxhi(i)*dzhiq(k)
         iout1(i, j, k) = iout1(i, j, k) + (utang1(i, j, k) - utang1(i, j, km))*emom*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5 !
       END DO

      !v
      !       il = block(n, 1)
      !       jl = MAX(block(n, 3) - myid*jmax, 1) + 1
      !       DO j = jl, ju
      !          DO i = il, iu
      !
      !             utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
      !             utang2Int = utang2(i, j, k)
      !             utangInt = max(umin, (utang1Int**2 + utang2Int**2))
      !             dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
      !             Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
      !             ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
      !             dummy = (utang2Int**2)*ctm
      !             bcmomflux = SIGN(dummy, utang2Int)
      !             iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
      !             eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
      !             iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k) !
      !          END DO
      !       END DO
      !
      ! !v top edge south
      !       j = MAX(block(n, 3) - myid*jmax, 1)
      !       DO i=il,iu
      !          utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
      !          utang2Int = utang2(i, j, k)
      !          utangInt = max(umin, (utang1Int**2 + utang2Int**2))
      !          dT = Tcell(i, j, k) - Twall
      !          Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri
      !          dummy = (utang2Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
      !          bcmomflux = SIGN(dummy, utang2Int)
      !          iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)*0.5
      !          eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
      !          iout2(i, j, k) = iout2(i, j, k) + ((utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k))*0.5 !
      !       END DO
      !
      ! !v top edge north
      !       j = MIN(block(n, 4) - myid*jmax, jmax) + 1
      !
      !       DO i = il, iu
      !          utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
      !          utang2Int = utang2(i, j, k)
      !          utangInt = max(umin, (utang1Int**2 + utang2Int**2))
      !          dT = Tcell(i, j - 1, k) - Twall
      !          Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri
      !          dummy = (utang2Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
      !          bcmomflux = SIGN(dummy, utang2Int)
      !          iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)*0.5
      !          eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
      !          iout2(i, j, k) = iout2(i, j, k) + ((utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k))*0.5 !
      !       END DO

      jl = block(n, 3)
      ju = block(n, 4)
      il = block(n, 1)
      iu = block(n, 2)

      if (jl>je+je*myid .or. ju < jb+je*myid) THEN

      else if (jl >= jb+je*myid .and. ju <= je+je*myid) then
        jl = block(n,3) - myid*je
        ju = block(n,4) - myid*je
        if (jl == ju) THEN
          do i = il,iu
            if(IIc(i,jl-1,k-1)==0 .and. IIc(i,jl-1,k)==1) then
              j = jl
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
              dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j-1, k)) + Tcell(i, j-1, k)*IIc(i, j-1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)
            else if (IIc(i,jl-1,k-1)==1 .or. IIc(i,jl-1,k) == 0 ) then
              j = jl
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = Tcell(i, j, k) - Twall
              dT = (Tcell(i, j, k) - Twall)*IIc(i,j,k)
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
            end if
            if(IIc(i,jl+1,k-1)==1 .and. IIc(i,jl+1,k)==1) then
              j = jl+1
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = Tcell(i, j-1, k) - Twall
              dT = (Tcell(i, j-1, k) - Twall)*IIc(i,j-1,k)
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
            else if (IIc(i,jl+1,k-1)==0 .and. IIc(i,jl+1,k)==1) then
            end if
          end do
        else if (jl /= ju) then
          do i = il,iu
            do j = jl+1, ju
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
              dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j-1, k)) + Tcell(i, j-1, k)*IIc(i, j-1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k) !
            end do
            if (IIc(i,jl-1,k) == 0 .or. IIc(i,jl-1,k-1) == 1) then
              j = jl
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = Tcell(i, j, k) - Twall
              dT = (Tcell(i, j, k) - Twall)*IIc(i,j,k)
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + ((utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k))*0.5
            else if (IIc(i,jl-1,k) == 1 .and. IIc(i,jl-1,k-1)== 0) then
              j = jl
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
              dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j-1, k)) + Tcell(i, j-1, k)*IIc(i, j-1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)
            end if
            if (IIc(i,ju+1,k) == 1 .and. IIc(i,ju+1,k-1) == 1) then
              j = ju+1
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = Tcell(i, j - 1, k) - Twall
              dT = (Tcell(i, j - 1, k) - Twall)*IIc(i,j-1,k)
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
            else
            end if
          end do
        end if
      else if (jl < jb+myid*je .and. ju > je+je*myid) then
        jl = jb
        ju = je
        if (jl == ju) THEN
          do i = il,iu
            if(IIc(i,jl-1,k-1)==0 .and. IIc(i,jl-1,k)==1) then
              j = jl
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
              dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j-1, k)) + Tcell(i, j-1, k)*IIc(i, j-1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)
            else if (IIc(i,jl-1,k-1)==1 .or. IIc(i,jl-1,k) == 0 ) then
              j = jl
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = Tcell(i, j, k) - Twall
              dT = (Tcell(i, j, k) - Twall)*IIc(i,j,k)
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
            end if
          end do
        else if (jl /= ju) then
          DO i =il,iu
            DO j = jl+1, ju
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
              dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j-1, k)) + Tcell(i, j-1, k)*IIc(i, j-1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k) !
            END DO
            if (IIc(i,jl-1,k-1)==0 .and. IIc(i,jl-1,k)==1) then
              j = jl
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
              dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j-1, k)) + Tcell(i, j-1, k)*IIc(i, j-1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k) !
            else if (IIc(i,jl-1,k-1)==1 .or. IIc(i,jl-1,k)==0) then
              j = jl
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = Tcell(i, j, k) - Twall
              dT = (Tcell(i, j, k) - Twall)*IIc(i,j,k)
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
            end if
          END DO
        end if
      else if (ju<=je+myid*je .and. jl<jb+myid*je) then
        jl = jb
        ju = block(n,4) - myid*je
        if (jl == ju) then
          do i = il,iu
            if(IIc(i,jl-1,k-1)==0 .and. IIc(i,jl-1,k)==1) then
              j = jl
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
              dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j-1, k)) + Tcell(i, j-1, k)*IIc(i, j-1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)
            else if (IIc(i,jl-1,k) == 0 ) then
              j = jl
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = Tcell(i, j, k) - Twall
              dT = (Tcell(i, j, k) - Twall)*IIc(i,j,k)
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
            end if
            if(IIc(i,jl+1,k-1)==1 .or. IIc(i,jl+1,k)==0) then
              j = jl+1
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = Tcell(i, j-1, k) - Twall
              dT = (Tcell(i, j-1, k) - Twall)*IIc(i,j-1,k)
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
            else if (IIc(i,jl+1,k-1)==0 .and. IIc(i,jl+1,k)==1) then
            end if
          end do
        else if (jl /= ju) then
          do i = il,iu
            do j = jl+1, ju
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
              dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j-1, k)) + Tcell(i, j-1, k)*IIc(i, j-1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k) !
            end do
            if (IIc(i,jl-1,k)==0) then
              j = jl
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = Tcell(i, j, k) - Twall
              dT = (Tcell(i, j, k) - Twall)*IIc(i,j,k)
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
            else if (IIc(i,jl-1,k)==1) then
              j = jl
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
              dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j-1, k)) + Tcell(i, j-1, k)*IIc(i, j-1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k) !
            end if
            if (IIc(i,ju+1,k-1)==1 .and. IIc(i,ju+1,k)==1) then
              j = ju+1
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = Tcell(i, j - 1, k) - Twall
              dT = (Tcell(i, j - 1, k) - Twall)*IIc(i,j-1,k)
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
            else if (IIc(i,ju+1,k)==0) THEN
            end if
          end do
        end if
      else if (ju>je+myid*je .and. jl>=jb+myid*je) then
        jl = block(n,3)-myid*je
        ju = je
        if (jl == ju) THEN
          j = jl
          do i = il,iu
            if(IIc(i,j-1,k-1) == 1 .or. IIc(i,j-1,k)==0) THEN
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = Tcell(i, j, k) - Twall
              dT = (Tcell(i, j, k) - Twall)*IIc(i,j,k)
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
            else if (IIc(i,j-1,k-1) == 0 .and. IIc(i,j-1,k)==1) THEN
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
              dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j-1, k)) + Tcell(i, j-1, k)*IIc(i, j-1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k) !
            end if
          end do
        else if (jl /= ju) then
          do i = il,iu
            do j = jl+1, ju
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
              dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j-1, k)) + Tcell(i, j-1, k)*IIc(i, j-1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k) !
            end do
            if (IIc(i,jl-1,k-1)==1 .or. IIc(i,jl-1,k)==0) then
              j = jl
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = Tcell(i, j, k)  - Twall
              dT = (Tcell(i, j, k)  - Twall)*IIc(i,j,k)
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
            else if (IIc(i,jl-1,k-1)==0 .and. IIc(i,jl-1,k)==1) then
              j = jl
              utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
              utang2Int = utang2(i, j, k)
              utangInt = max(umin, (utang1Int**2 + utang2Int**2))
              !dT = Tcell(i, j, k)  - Twall
              dT = (Tcell(i, j, k)  - Twall)*IIc(i,j,k)
              Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
              ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
              !dummy = (utang2Int**2)*ctm
              dummy = abs(utang2Int)*sqrt(utangInt)*ctm
              bcmomflux = SIGN(dummy, utang2Int)
              iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
              eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)
              iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
            end if
          end do
        end if
      end if
      !end if

   ! CASE (61) !wall in xy -> wf in z (=horizontal), bottom wall
   !    k = block(n, 5) - 1 !block location
   !    !if (.not.(k.lt.0)) then
   !    if (k < kb) return
   !    kp = k + 1 !shear velocity location
   !    il = block(n, 1) + 1
   !    iu = block(n, 2)
   !    jl = MAX(block(n, 3) - myid*jmax, 1)
   !    ju = MIN(block(n, 4) - myid*jmax, jmax)
   !
   !    delta = 0.5*dzf(k)
   !    logdz = LOG(delta/z0)
   !    logdzh = LOG(delta/z0h)
   !    logzh = LOG(z0/z0h)
   !    sqdz = SQRT(delta/z0)
   !
   !    DO j = jl, ju
   !      DO i = il, iu
   !        ! utang1Int = utang1(i, j, k)
   !        ! utang2Int = (utang2(i, j, k) + utang2(i - 1, j, k) + utang2(i, j + 1, k) + utang2(i - 1, j + 1, k))*0.25
   !        ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !        ! dT = ((Tcell(i, j, k) + Tcell(i - 1, j, k)) - (Twall + Twall))*0.5
   !        ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !        ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !        ! dummy = (utang1Int**2)*ctm
   !        ! bcmomflux = SIGN(dummy, utang1Int)
   !        ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !        emop = ( dzf(kp) * ( ekm(i,j,k )*dxf(i-1) + ekm(i-1,j,k )*dxf(i) )  + &
   !        dzf(k)  * ( ekm(i,j,kp)*dxf(i-1) + ekm(i-1,j,kp)*dxf(i) ) )*dxhi(i) * dzhiq(kp)
   !        iout1(i, j, k) = iout1(i, j, k) - (utang1(i, j, kp) - utang1(i, j, k))*emop*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k) !
   !      END DO
   !    END DO
   !
   !    !u top edge west
   !    i = block(n, 1)
   !    DO j = jl, ju
   !      ! utang1Int = utang1(i, j, k)
   !      ! utang2Int = (utang2(i, j, k) + utang2(i - 1, j, k) + utang2(i, j + 1, k) + utang2(i - 1, j + 1, k))*0.25
   !      ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !      ! dT = Tcell(i, j, k) - Twall
   !      ! Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri
   !      ! dummy = (utang1Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !      ! bcmomflux = SIGN(dummy, utang1Int)
   !      ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)*0.5
   !      emop = ( dzf(kp) * ( ekm(i,j,k )*dxf(i-1) + ekm(i-1,j,k )*dxf(i) )  + &
   !      dzf(k)  * ( ekm(i,j,kp)*dxf(i-1) + ekm(i-1,j,kp)*dxf(i) ) )*dxhi(i) * dzhiq(kp)
   !      iout1(i, j, k) = iout1(i, j, k) - (utang1(i, j, kp) - utang1(i, j, k))*emop*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5 !
   !    END DO
   !
   !    !u top edge east
   !    i = block(n, 2)
   !    DO j = jl, ju
   !      ! utang1Int = utang1(i, j, k)
   !      ! utang2Int = (utang2(i, j, k) + utang2(i - 1, j, k) + utang2(i, j + 1, k) + utang2(i - 1, j + 1, k))*0.25
   !      ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !      ! dT = Tcell(i - 1, j, k) - Twall
   !      ! Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri
   !      ! dummy = (utang1Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !      ! bcmomflux = SIGN(dummy, utang1Int)
   !      ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)*0.5
   !      emop = ( dzf(kp) * ( ekm(i,j,k )*dxf(i-1) + ekm(i-1,j,k )*dxf(i) )  + &
   !      dzf(k)  * ( ekm(i,j,kp)*dxf(i-1) + ekm(i-1,j,kp)*dxf(i) ) )*dxhi(i) * dzhiq(kp)
   !      iout1(i, j, k) = iout1(i, j, k) - (utang1(i, j, kp) - utang1(i, j, k))*emop*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5 !
   !    END DO
   !
   !    ! !v
   !    !    il = block(n, 1)
   !    !    jl = MAX(block(n, 3) - myid*jmax, 1) + 1
   !    !    DO j = jl, ju
   !    !       DO i = il, iu
   !    !          ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !    !          ! utang2Int = utang2(i, j, k)
   !    !          ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !    !          ! dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
   !    !          ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !    !          ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !    !          ! dummy = (utang2Int**2)*ctm
   !    !          ! bcmomflux = SIGN(dummy, utang2Int)
   !    !          ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !    !          eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !    !                   dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !    !          iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k) !
   !    !       END DO
   !    !    END DO
   !    !
   !    ! !v top edge south
   !    !    j = MAX(block(n, 3) - myid*jmax, 1)
   !    !    DO i=il,iu
   !    !       ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !    !       ! utang2Int = utang2(i, j, k)
   !    !       ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !    !       ! dT = Tcell(i, j, k) - Twall
   !    !       ! Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri
   !    !       ! dummy = (utang2Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !    !       ! bcmomflux = SIGN(dummy, utang2Int)
   !    !       ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)*0.5
   !    !       eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !    !                dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !    !       iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5 !
   !    !    END DO
   !    !
   !    ! !v top edge north
   !    !    j = MIN(block(n, 4) - myid*jmax, jmax) + 1
   !    !
   !    !    DO i = il, iu
   !    !       ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !    !       ! utang2Int = utang2(i, j, k)
   !    !       ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !    !       ! dT = Tcell(i, j - 1, k) - Twall
   !    !       ! Ribl0 = grav*delta*dT/(Twall*utangInt) !Eq. 6, guess initial Ri
   !    !       ! dummy = (utang2Int**2)*unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !    !       ! bcmomflux = SIGN(dummy, utang2Int)
   !    !       ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)*0.5
   !    !       eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !    !                dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !    !       iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5 !
   !    !    END DO
   !
   !    jl = block(n, 3)
   !    ju = block(n, 4)
   !    il = block(n, 1)
   !    iu = block(n, 2)
   !
   !
   !    if (jl>je+je*myid .or. ju < jb+je*myid) THEN
   !
   !    else if (jl >= jb+je*myid .and. ju <= je+je*myid) then
   !      jl = block(n,3) - myid*je
   !      ju = block(n,4) - myid*je
   !      if (jl == ju) THEN
   !        do i = il,iu
   !          if(IIc(i,jl-1,k+1)==0 .and. IIc(i,jl-1,k)==1) then
   !            j = jl
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
   !            ! dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j-1, k)) + Tcell(i, j-1, k)*IIc(i, j-1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)
   !          else if (IIc(i,jl-1,k+1)==1 .or. IIc(i,jl-1,k) == 0 ) then
   !            j = jl
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = Tcell(i, j, k) - Twall
   !            ! dT = (Tcell(i, j, k) - Twall)*IIc(i,j,k)
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
   !          end if
   !          if(IIc(i,jl+1,k+1)==1 .and. IIc(i,jl+1,k)==1) then
   !            j = jl+1
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = Tcell(i, j-1, k) - Twall
   !            ! dT = (Tcell(i, j-1, k) - Twall)*IIc(i,j-1,k)
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  ) + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
   !          else if (IIc(i,jl+1,k+1)==0 .and. IIc(i,jl+1,k)==1) then
   !          end if
   !        end do
   !      else if (jl /= ju) then
   !        do i = il,iu
   !          do j = jl+1, ju
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
   !            ! dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j-1, k)) + Tcell(i, j-1, k)*IIc(i, j-1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)
   !          end do
   !          if (IIc(i,jl-1,k) == 0 .or. IIc(i,jl-1,k+1) == 1) then
   !            j = jl
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = Tcell(i, j, k) - Twall
   !            ! dT = (Tcell(i, j, k) - Twall)*IIc(i,j,k)
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
   !          else if (IIc(i,jl-1,k) == 1 .and. IIc(i,jl-1,k+1)== 0) then
   !            j = jl
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
   !            ! dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j-1, k)) + Tcell(i, j-1, k)*IIc(i, j-1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)
   !          end if
   !          if (IIc(i,ju+1,k) == 1 .and. IIc(i,ju+1,k+1) == 1) then
   !            ! j = ju+1
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = Tcell(i, j - 1, k) - Twall
   !            ! dT = (Tcell(i, j - 1, k) - Twall)*IIc(i,j-1,k)
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
   !          else
   !          end if
   !        end do
   !      end if
   !    else if (jl < jb+myid*je .and. ju > je+je*myid) then
   !      jl = jb
   !      ju = je
   !      if (jl == ju) THEN
   !        do i = il,iu
   !          if(IIc(i,jl-1,k+1)==0 .and. IIc(i,jl-1,k)==1) then
   !            j = jl
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
   !            ! dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j-1, k)) + Tcell(i, j-1, k)*IIc(i, j-1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)
   !          else if (IIc(i,jl-1,k+1)==1 .or. IIc(i,jl-1,k) == 0 ) then
   !            j = jl
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = Tcell(i, j, k) - Twall
   !            ! dT = (Tcell(i, j, k) - Twall)*IIc(i,j,k)
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
   !          end if
   !        end do
   !      else if (jl /= ju) then
   !        DO i =il,iu
   !          DO j = jl+1, ju
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
   !            ! dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j-1, k)) + Tcell(i, j-1, k)*IIc(i, j-1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)
   !          END DO
   !          if (IIc(i,jl-1,k+1)==0 .and. IIc(i,jl-1,k)==1) then
   !            j = jl
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
   !            ! dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j-1, k)) + Tcell(i, j-1, k)*IIc(i, j-1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)
   !          else if (IIc(i,jl-1,k+1)==1 .or. IIc(i,jl-1,k)==0) then
   !            j = jl
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = Tcell(i, j, k) - Twall
   !            ! dT = (Tcell(i, j, k) - Twall)*IIc(i,j,k)
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
   !          end if
   !        END DO
   !      end if
   !    else if (ju<=je+myid*je .and. jl<jb+myid*je) then
   !      jl = jb
   !      ju = block(n,4) - myid*je
   !      if (jl == ju) then
   !        do i = il,iu
   !          if(IIc(i,jl-1,k-1)==0 .and. IIc(i,jl-1,k)==1) then
   !            j = jl
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
   !            ! dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j-1, k)) + Tcell(i, j-1, k)*IIc(i, j-1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)
   !          else if (IIc(i,jl-1,k) == 0 ) then
   !            j = jl
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = Tcell(i, j, k) - Twall
   !            ! dT = (Tcell(i, j, k) - Twall)*IIc(i,j,k)
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
   !          end if
   !          if(IIc(i,jl+1,k+1)==1 .or. IIc(i,jl+1,k)==0) then
   !            j = jl+1
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = Tcell(i, j-1, k) - Twall
   !            ! dT = (Tcell(i, j-1, k) - Twall)*IIc(i,j-1,k)
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
   !          else if (IIc(i,jl+1,k-1)==0 .and. IIc(i,jl+1,k)==1) then
   !          end if
   !        end do
   !      else if (jl /= ju) then
   !        do i = il,iu
   !          do j = jl+1, ju
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
   !            ! dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j-1, k)) + Tcell(i, j-1, k)*IIc(i, j-1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)
   !          end do
   !          if (IIc(i,jl-1,k)==0) then
   !            j = jl
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = Tcell(i, j, k) - Twall
   !            ! dT = (Tcell(i, j, k) - Twall)*IIc(i,j,k)
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
   !          else if (IIc(i,jl-1,k)==1) then
   !            j = jl
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
   !            ! dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j-1, k)) + Tcell(i, j-1, k)*IIc(i, j-1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)
   !          end if
   !          if (IIc(i,ju+1,k+1)==1 .and. IIc(i,ju+1,k)==1) then
   !            j = ju+1
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = Tcell(i, j - 1, k) - Twall
   !            ! dT = (Tcell(i, j - 1, k) - Twall)*IIc(i,j-1,k)
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
   !          else if (IIc(i,ju+1,k)==0) THEN
   !          end if
   !        end do
   !      end if
   !    else if (ju>je+myid*je .and. jl>=jb+myid*je) then
   !      jl = block(n,3)-myid*je
   !      ju = je
   !      if (jl == ju) THEN
   !        j = jl
   !        do i = il,iu
   !          if(IIc(i,j-1,k+1) == 1 .or. IIc(i,j-1,k)==0) THEN
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = Tcell(i, j, k) - Twall
   !            ! dT = (Tcell(i, j, k) - Twall)*IIc(i,j,k)
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
   !          else if (IIc(i,j-1,k+1) == 0 .and. IIc(i,j-1,k)==1) THEN
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
   !            ! dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j-1, k)) + Tcell(i, j-1, k)*IIc(i, j-1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)
   !          end if
   !        end do
   !      else if (jl /= ju) then
   !        do i = il,iu
   !          do j = jl+1, ju
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
   !            ! dT = (Tcell(i, j, k)*IIc(i, j, k)*(2-IIc(i, j-1, k)) + Tcell(i, j-1, k)*IIc(i, j-1, k)*(2-IIc(i, j, k)) - (Twall + Twall))*0.5
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)
   !          end do
   !          if (IIc(i,jl-1,k+1)==1 .or. IIc(i,jl-1,k)==0) then
   !            j = jl
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = Tcell(i, j, k)  - Twall
   !            ! dT = (Tcell(i, j, k)  - Twall)*IIc(i,j,k)
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
   !          else if (IIc(i,jl-1,k+1)==0 .and. IIc(i,jl-1,k)==1) then
   !            j = jl
   !            ! utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
   !            ! utang2Int = utang2(i, j, k)
   !            ! utangInt = max(umin, (utang1Int**2 + utang2Int**2))
   !            ! !dT = Tcell(i, j, k)  - Twall
   !            ! dT = (Tcell(i, j, k)  - Twall)*IIc(i,j,k)
   !            ! Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
   !            ! ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
   !            ! dummy = (utang2Int**2)*ctm
   !            ! bcmomflux = SIGN(dummy, utang2Int)
   !            ! iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
   !            eomp = ( dzf(kp) * ( ekm(i,j,k ) + ekm(i,jm,k )  )  + &
   !            dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * dzhiq(kp)
   !            iout2(i, j, k) = iout2(i, j, k) - (utang2(i, j, kp) - utang2(i, j, k))*eomp*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k)*0.5
   !          end if
   !        end do
   !      end if
   !    end if
   ! !end if

!!!!!!!!!!!!!!!SPECIAL CASES FOR THE SURFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!can actually be made redundant and just be replaced by standard horizontal case (doesn't really matter though)
   CASE (91) !surface momentum flux

      k = kb !
      km = k - 1 !
      il = ib
      iu = ie
      jl = jb
      ju = je


      delta = 0.5*dzf(k) !might need attention on streched grids! as well as the dzfi when updating up
            logdz = LOG(delta/z0)
            logdzh = LOG(delta/z0h)
            logzh = LOG(z0/z0h)
            sqdz = SQRT(delta/z0)

      DO j = jl, ju !u component
         DO i = il, iu
            utang1Int = utang1(i, j, k)
            utang2Int = (utang2(i, j, k) + utang2(i - 1, j, k) + utang2(i, j + 1, k) + utang2(i - 1, j + 1, k))*0.25
            utangInt = max(umin, (utang1Int**2 + utang2Int**2))
            dT = ((Tcell(i, j, k) + Tcell(i - 1, j, k)) - (Twall + Twall))*0.5
            Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
            ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
            !dummy = (utang1Int**2)*ctm
            dummy = abs(utang1Int)*sqrt(utangInt)*ctm
            bcmomflux = SIGN(dummy, utang1Int) !bcmomflux=u_star^2

            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)

            emom = (dzf(km)*(ekm(i, j, k)*dxf(i - 1) + ekm(i - 1, j, k)*dxf(i)) + & ! dx is non-equidistant
                    dzf(k)*(ekm(i, j, km)*dxf(i - 1) + ekm(i - 1, j, km)*dxf(i)))*dxhi(i)*dzhiq(k)

            iout1(i, j, k) = iout1(i, j, k) + (utang1(i, j, k) - utang1(i, j, km))*emom*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k) !

         END DO
      END DO

      DO j = jl, ju
         DO i = il, iu

            utang1Int = (utang1(i, j, k) + utang1(i, j - 1, k) + utang1(i + 1, j - 1, k) + utang1(i + 1, j, k))*0.25
            utang2Int = utang2(i, j, k)
            utangInt = max(umin, (utang1Int**2 + utang2Int**2))
            dT = ((Tcell(i, j, k) + Tcell(i, j - 1, k)) - (Twall + Twall))*0.5
            Ribl0 = grav*delta*dT*2/((Twall + Twall)*utangInt) !Eq. 6, guess initial Ri
            ctm = unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !save result and update field
            !dummy = (utang2Int**2)*ctm !save result and update field
            dummy = abs(utang2Int)*sqrt(utangInt)*ctm
            bcmomflux = SIGN(dummy, utang2Int)
            iomomflux(i, j, k) = iomomflux(i, j, k) + bcmomflux*dzfi(k)
            eomm = (dzf(km)*(ekm(i, j, k) + ekm(i, j - 1, k)) + dzf(k)*(ekm(i, j, km) + ekm(i, j - 1, km)))*dzhiq(k)

            iout2(i, j, k) = iout2(i, j, k) + (utang2(i, j, k) - utang2(i, j, km))*eomm*dzhi(k)*dzfi(k) - bcmomflux*dzfi(k) !
         END DO
      END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   CASE (92) !surface temperature flux
      k = kb !block location
      ku = k !shear velocity location
      il = ib
      iu = ie
      jl = jb
      ju = je

      delta = dzf(k)*0.5
            logdz = LOG(delta/z0)
            logdzh = LOG(delta/z0h)
            logzh = LOG(z0/z0h)
            sqdz = SQRT(delta/z0)

      DO j = jl, ju
         DO i = il, iu

            utang1Int = (utang1(i, j, ku) + utang1(i + 1, j, ku))*0.5
            utang2Int = (utang2(i, j, ku) + utang2(i, j + 1, ku))*0.5
            utangInt = max(umin, (utang1Int**2 + utang2Int**2))
            dT = (Tcell(i, j, ku) - Twall)

            Ribl0 = grav*delta*dT/(Twall*utangInt) !

            call unoh(bcTflux, cth, logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2)
            obcTfluxA = obcTfluxA + bcTflux
            iotflux(i, j, k) = iotflux(i, j, k) + bcTflux*dzfi(k)

            iot(i, j, ku) = iot(i, j, ku) + &
                            0.5*(dzf(k - 1)*ekh(i, j, k) + dzf(k)*ekh(i, j, k - 1))* & ! zero flux
                            (Tcell(i, j, k) - Tcell(i, j, k - 1))*dzh2i(ku)*dzfi(ku) &
                            - bcTflux*dzfi(k)

         END DO
      END DO

   END SELECT
END SUBROUTINE wfuno

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!for scalar
!FUNCTION unoh(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !for heat, the bit that does not change no matter what wall
SUBROUTINE unoh(otf, octh, logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !for heat, the bit that does not change no matter what wall
   use modglobal, only : prandtlturb
!flux in Km/s
   IMPLICIT NONE
   REAL, INTENT(in) :: logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2
   REAL, INTENT(out) :: octh, otf
   REAL :: Ribl1, Fm, Fh, cm, ch, M, dTrough, cth
   REAL, PARAMETER :: b1 = 9.4 !parameters from Uno1995
   REAL, PARAMETER :: b2 = 4.7
   REAL, PARAMETER :: dm = 7.4
   REAL, PARAMETER :: dh = 5.3

   octh = 0.
   otf = 0.
   IF (Ribl0 > 0.) THEN !0.25 approx critical for bulk Richardson number  => stable
      Fm = 1./(1. + b2*Ribl0)**2 !Eq. 4
      Fh = Fm !Eq. 4
   ELSE ! => unstable
      cm = (dm*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
      ch = (dh*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
      Fm = 1. - (b1*Ribl0)/(1. + cm*SQRT(ABS(Ribl0))) !Eq. 3
      Fh = 1. - (b1*Ribl0)/(1. + ch*SQRT(ABS(Ribl0))) !Eq. 3
   END IF

   M = prandtlturb*logdz*SQRT(Fm)/Fh !Eq. 14

   Ribl1 = Ribl0 - Ribl0*prandtlturb*logzh/(prandtlturb*logzh + M) !Eq. 17

   !interate to get new Richardson number
   IF (Ribl1 > 0.) THEN !0.25 approx critical for bulk Richardson number  => stable
      Fm = 1./(1. + b2*Ribl1)**2 !Eq. 4
      Fh = Fm !Eq. 4
   ELSE ! => unstable
      cm = (dm*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
      ch = (dh*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
      Fm = 1. - (b1*Ribl1)/(1. + cm*SQRT(ABS(Ribl1))) !Eq. 3
      Fh = 1. - (b1*Ribl1)/(1. + ch*SQRT(ABS(Ribl1))) !Eq. 3
   END IF

   ! ! Uno (2)
   ! M = prandtlturb*logdz*SQRT(Fm)/Fh !Eq. 14
   ! dTrough = dT*1./(prandtlturb*logzh/M + 1.) !Eq. 13a
   ! octh = SQRT(utangInt)*fkar2/(logdz*logdz)*Fh/prandtlturb !Eq. 8
   ! otf = octh*dTrough !Eq. 2, Eq. 8

   ! Uno (8)
   octh = SQRT(utangInt)*fkar2/(logdz*logdzh)*Fh/prandtlturb !Eq. 8
   otf = octh*dT !Eq. 2, Eq. 8

END SUBROUTINE unoh
!END FUNCTION unoh

!!!!!!!!!!!!!
!for momentum
REAL FUNCTION unom(logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2) !for momentum, this bit is not depended on orientation etc
   use modglobal, only : prandtlturb
!momentum flux in m2/s2
!dT,utang and logdzh are unused and could be removed
   IMPLICIT NONE
   REAL, INTENT(in) :: logdz, logdzh, logzh, sqdz, utangInt, dT, Ribl0, fkar2
   REAL :: Ribl1, Fm, Fh, cm, ch, Ctm, M
   REAL, PARAMETER :: b1 = 9.4 !parameters from Uno1995
   REAL, PARAMETER :: b2 = 4.7
   REAL, PARAMETER :: dm = 7.4
   REAL, PARAMETER :: dh = 5.3
   !REAL, PARAMETER :: prandtlmol = 0.71

   IF (Ribl0 > 0.) THEN !0.25 approx critical for bulk Richardson number  => stable
      Fm = 1./(1. + b2*Ribl0)**2 !Eq. 4
      Fh = Fm !Eq. 4
   ELSE ! => unstable
      cm = (dm*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
      ch = (dh*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
      Fm = 1. - (b1*Ribl0)/(1. + cm*SQRT(ABS(Ribl0))) !Eq. 3
      Fh = 1. - (b1*Ribl0)/(1. + ch*SQRT(ABS(Ribl0))) !Eq. 3
   END IF

   M = prandtlturb*logdz*SQRT(Fm)/Fh !Eq. 14

   Ribl1 = Ribl0 - Ribl0*prandtlturb*logzh/(prandtlturb*logzh + M) !Eq. 17

   !interate to get new Richardson number
   IF (Ribl1 > 0.) THEN !0.25 approx critical for bulk Richardson number  => stable
      Fm = 1./(1. + b2*Ribl1)**2 !Eq. 4
   ELSE ! => unstable
      cm = (dm*fkar2)/(logdz**2)*b1*sqdz !Eq. 5
      Fm = 1. - (b1*Ribl1)/(1. + cm*SQRT(ABS(Ribl1))) !Eq. 3
   END IF

   Ctm = fkar2/(logdz**2)*Fm !Eq. 7
   unom = Ctm !Eq. 2, Eq. 8
END FUNCTION unom
