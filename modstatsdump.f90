!> \file modstatsdump.f90
!!  Dumps statistics of various fields
!>
!!  \author Tom Grylls, ICL May 25 2016
!
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
module modstatsdump

  use modglobal, only : dt,lydump,lytdump,ltkedump,lxydump,lxytdump,ltdump,ifoutput !,nstat
  use modmpi, only : myid
  implicit none
  private
  PUBLIC :: initstatsdump,statsdump,exitstatsdump
  save

  !NetCDF variables
  integer :: ncidy,ncidyt,ncidtke,ncidxy,ncidslice,ncidxyt,nrecy=0,nrecyt=0,nrectke=0,nrecxy=0,&
             nrecslice=0,nrecxyt=0,nstatyt=14,nstaty=14,nstattke=8,nstatxy=14,nstatslice=8,&
             nstatxyt=19,ncidt,nrect=0,nstatt=15
  character(80) :: yname = 'ydump.xxx.nc'
  character(80) :: ytname = 'ytdump.xxx.nc'
  character(80) :: tkename = 'tkedump.xxx.nc'
  character(80) :: xyname = 'xydump.xxx.nc'
  character(80) :: xytname = 'xytdump.xxx.nc'
  character(80) :: tname = 'tdump.xxx.xxx.nc'
  character(80) :: slicename = 'slicedump.xxx.xxx.nc'
  character(80),dimension(1,4) :: tncstaty
  character(80),dimension(1,4) :: tncstatyt
  character(80),dimension(1,4) :: tncstattke
  character(80),dimension(1,4) :: tncstatxy
  character(80),dimension(1,4) :: tncstatslice
  character(80),dimension(1,4) :: tncstatxyt
  character(80),dimension(1,4) :: tncstatt

  integer :: klow,khigh,i,j,k
  real    :: tsamplep,tstatsdumpp,tsample,tstatsdump

contains

  !--------------------------
  !> Initializing statsdump. Read out the namelist, initializing the variables
  !-------------------------

  subroutine initstatsdump
    use modmpi,   only : my_real,mpierr,comm3d,mpi_logical,mpi_integer,mpi_character,cmyid
    use modglobal,only : imax,jmax,kmax,cexpnr,ifnamopt,fname_options,kb,ke,ladaptive,btime,&
                         nsv,lslicedump,lxytdump
    use modstat_nc,only: open_nc, define_nc,ncinfo,writestat_dims_nc
    use modfields,only : ncstaty,ncstatyt,ncstattke,ncstatxy,ncstatslice,ncstatxyt,ncstatt
    implicit none
    integer :: ierr

    namelist/NAMSTATSDUMP/ &
         lydump,tsample,klow,khigh,tstatsdump,lytdump,ltkedump,lxydump,lxytdump,ltdump

    allocate(ncstaty(nstaty,4))
    allocate(ncstatyt(nstatyt,4))
    allocate(ncstattke(nstattke,4))
    allocate(ncstatxy(nstatxy,4))
    allocate(ncstatslice(nstatslice,4))
    allocate(ncstatxyt(nstatxyt,4))
    allocate(ncstatt(nstatt,4))

    klow=kb
    khigh=ke

    if(myid==0)then
       open(ifnamopt,file=fname_options,status='old',iostat=ierr)
       read (ifnamopt,NAMSTATSDUMP,iostat=ierr)
       if (ierr > 0) then
          print *, 'Problem in namoptions NAMSTATSDUMP'
          print *, 'iostat error: ', ierr
          stop 'ERROR: Problem in namoptions NAMSTATSDUMP'
       endif
       write(6 ,NAMSTATSDUMP)
       close(ifnamopt)
    end if

    call MPI_BCAST(klow        ,1,MPI_INTEGER,0,comm3d,ierr) !have to do this? just want nc for first CPU
    call MPI_BCAST(khigh       ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(nstatt      ,1,MPI_INTEGER,0,comm3d,ierr)
!    call MPI_BCAST(nstaty      ,1,MPI_INTEGER,0,comm3d,ierr)
    call MPI_BCAST(ncstatyt    ,80,MPI_CHARACTER,0,comm3d,mpierr)
    call MPI_BCAST(ncstaty     ,80,MPI_CHARACTER,0,comm3d,mpierr)
    call MPI_BCAST(ncstattke   ,80,MPI_CHARACTER,0,comm3d,mpierr)
    call MPI_BCAST(ncstatxy    ,80,MPI_CHARACTER,0,comm3d,mpierr)
    call MPI_BCAST(ncstatxyt   ,80,MPI_CHARACTER,0,comm3d,mpierr)
    call MPI_BCAST(ncstatt     ,80,MPI_CHARACTER,0,comm3d,mpierr)
    call MPI_BCAST(ltdump      ,1,MPI_LOGICAL,0,comm3d,ierr)

    !> Generate y-averaged NetCDF: ydump.xxx.nc  
    if(lydump) then

      yname(7:9) = cexpnr
      call ncinfo(tncstaty(1,:),'time'      ,'Time'                         ,'s'      ,'time')
      call ncinfo(ncstaty( 1,:),'uy'        ,'Streamwise velocity'          ,'m/s'    ,'m0tt')
      call ncinfo(ncstaty( 2,:),'vy'        ,'Spanwise velocity'            ,'m/s'    ,'t0tt')
      call ncinfo(ncstaty( 3,:),'wy'        ,'Vertical velocity'            ,'m/s'    ,'t0mt')
      call ncinfo(ncstaty( 4,:),'thly'      ,'Temperature'                  ,'K'      ,'t0tt')
      call ncinfo(ncstaty( 5,:),'qty'       ,'Moisture'                     ,'kg/kg'  ,'t0tt')
      call ncinfo(ncstaty( 6,:),'sca1y'     ,'Scalar field 1'               ,'kg/m^3' ,'t0tt')
      call ncinfo(ncstaty( 7,:),'sca2y'     ,'Scalar field 2'               ,'kg/m^3' ,'t0tt')
      call ncinfo(ncstaty( 8,:),'sca3y'     ,'Scalar field 3'               ,'kg/m^3' ,'t0tt')
      call ncinfo(ncstaty( 9,:),'upwpy'     ,'Turbulent mom. flux'          ,'m^2/s^2','m0mt')
      call ncinfo(ncstaty(10,:),'wpthlpy'   ,'Turbulent heat flux'          ,'K m/s'  ,'t0mt')
      call ncinfo(ncstaty(11,:),'usgsy'     ,'SGS mom. flux'                ,'m^2/s^2','m0mt') 
      call ncinfo(ncstaty(12,:),'thlsgsy'   ,'SGS heat flux'                ,'K m/s'  ,'t0mt')
      call ncinfo(ncstaty(13,:),'uwyik'     ,'Advective mom. flux'          ,'m^2/s^2','m0mt') 
      call ncinfo(ncstaty(14,:),'wthlyk'    ,'Advective heat flux'          ,'K m/s'  ,'t0mt')

      if (myid==0) then
        call open_nc(yname, ncidy, nrecy, n1=imax, n3=khigh-klow+1)
        if (nrecy==0) then
          call define_nc( ncidy, 1, tncstaty)
          call writestat_dims_nc(ncidy)
        end if
      call define_nc( ncidy, nstaty, ncstaty)
      endif !myid==0
    endif

    !> Generate time and y averaged NetCDF: ytdump.xxx.nc
    if (lytdump) then
    
      ytname(8:10) = cexpnr
      call ncinfo(tncstatyt(1,:),'time'       ,'Sampling time'             ,'s'       ,'time')
      call ncinfo(ncstatyt( 1,:),'uyt'        ,'Streamwise velocity'       ,'m/s'     ,'m0tt')
      call ncinfo(ncstatyt( 2,:),'vyt'        ,'Spanwise velocity'         ,'m/s'     ,'t0tt')
      call ncinfo(ncstatyt( 3,:),'wyt'        ,'Vertical velocity'         ,'m/s'     ,'t0mt')
      call ncinfo(ncstatyt( 4,:),'thlyt'      ,'Temperature'               ,'K'       ,'t0tt')
      call ncinfo(ncstatyt( 5,:),'qtyt'       ,'Moisture'                  ,'kg/kg'   ,'t0tt')
      call ncinfo(ncstatyt( 6,:),'sca1yt'     ,'Scalar field 1'            ,'kg/m^3'  ,'t0tt')
      call ncinfo(ncstatyt( 7,:),'sca2yt'     ,'Scalar field 2'            ,'kg/m^3'  ,'t0tt')
      call ncinfo(ncstatyt( 8,:),'sca3yt'     ,'Scalar field 3'            ,'kg/m^3'  ,'t0tt')
      call ncinfo(ncstatyt( 9,:),'upwpyt'     ,'Turbulent mom. flux'       ,'m^2/s^2' ,'m0mt')
      call ncinfo(ncstatyt( 10,:),'wpthlpyt'  ,'Turbulent heat flux'       ,'K m/s'   ,'t0mt')
      call ncinfo(ncstatyt( 11,:),'uwyt'      ,'Kinematic mom. flux'       ,'m^2/s^2' ,'m0mt')
      call ncinfo(ncstatyt( 12,:),'wthlyt'    ,'Kinematic heat flux'       ,'K m/s'   ,'t0mt')
      call ncinfo(ncstatyt( 13,:),'usgsyt'    ,'SGS mom. flux'             ,'m^2/s^2' ,'m0mt')
      call ncinfo(ncstatyt( 14,:),'thlsgsyt'  ,'SGS heat flux'             ,'K m/s'   ,'t0mt')

      if (myid==0) then
        call open_nc(ytname, ncidyt, nrecyt, n1=imax, n3=khigh-klow+1)
        if (nrecyt==0) then
          call define_nc( ncidyt, 1, tncstatyt)
          call writestat_dims_nc(ncidyt)
        end if
        call define_nc( ncidyt, nstatyt, ncstatyt)
      endif !myid==0
    endif

    !> Generate time, y and x averaged NetCDF: xydump.xxx.nc
    if (lxydump) then
    
      xyname(8:10) = cexpnr
      call ncinfo(tncstatxy(1,:),'time'    ,'Time'                        ,'s'      ,'time')
      call ncinfo(ncstatxy( 1,:),'uxy'     ,'Streamwise velocity'         ,'m/s'    ,'tt'  )
      call ncinfo(ncstatxy( 2,:),'vxy'     ,'Spanwise velocity'           ,'m/s'    ,'tt'  )
      call ncinfo(ncstatxy( 3,:),'wxy'     ,'Vertical velocity'           ,'m/s'    ,'mt'  )
      call ncinfo(ncstatxy( 4,:),'thlxy'   ,'Temperature'                 ,'K'      ,'tt'  )
      call ncinfo(ncstatxy( 5,:),'qtxy'    ,'Moisture'                    ,'kg/kg'  ,'tt'  )
      call ncinfo(ncstatxy( 6,:),'upwpxy'  ,'Mom. flux'                   ,'m^2/s^2','mt'  )
      call ncinfo(ncstatxy( 7,:),'wpthlpxy','Heat flux'                   ,'Km/s'   ,'mt'  )
      call ncinfo(ncstatxy( 8,:),'vpwpxy'  ,'Mom. flux'                   ,'Km/s'   ,'mt'  )
      call ncinfo(ncstatxy( 9,:),'usgsxy'  ,'SGS mom. flux'               ,'m^2/s^2','mt'  )
      call ncinfo(ncstatxy(10,:),'thlsgsxy','SGS heat flux'               ,'Km/s'   ,'mt'  )
      call ncinfo(ncstatxy(11,:),'vsgsxy'  ,'SGS mom. flux'               ,'m^2/s^2','mt'  )
      call ncinfo(ncstatxy(12,:),'uwxyik'  ,'Advective mom. flux'         ,'m^2/s^2','mt') 
      call ncinfo(ncstatxy(13,:),'wthlxy'  ,'Advective heat flux'         ,'K m/s'  ,'mt')
      call ncinfo(ncstatxy(14,:),'vwxy'    ,'Advective mom. flux'         ,'m^2/s^2','mt')
      if (myid==0) then      
        call open_nc(xyname, ncidxy, nrecxy, n3=khigh-klow+1)
        if (nrecxy==0) then
          call define_nc( ncidxy, 1, tncstatxy)
          call writestat_dims_nc(ncidxy)
        end if
        call define_nc( ncidxy, nstatxy, ncstatxy)
      end if
    end if

    !> Generate time, y and x averaged NetCDF: xytdump.xxx.nc
    if (lxytdump) then
    
      xytname(9:11) = cexpnr
      call ncinfo(tncstatxyt(1,:),'time'      ,'Time'                        ,'s'      ,'time')
      call ncinfo(ncstatxyt( 1,:),'uxyt'       ,'Streamwise velocity'         ,'m/s'    ,'tt'  )
      call ncinfo(ncstatxyt( 2,:),'vxyt'       ,'Spanwise velocity'           ,'m/s'    ,'tt'  )
      call ncinfo(ncstatxyt( 3,:),'wxyt'       ,'Vertical velocity'           ,'m/s'    ,'mt'  )
      call ncinfo(ncstatxyt( 4,:),'thlxyt'     ,'Temperature'                 ,'K'      ,'tt'  )
      call ncinfo(ncstatxyt( 5,:),'qtxyt'      ,'Moisture'                    ,'kg/kg'  ,'tt'  )
      call ncinfo(ncstatxyt( 6,:),'upwpxyt'    ,'Turbulent mom. flux'         ,'m^2/s^2','mt'  )
      call ncinfo(ncstatxyt( 7,:),'wpthlpxyt'  ,'Turbulent heat flux'         ,'K m/s'  ,'mt'  )
      call ncinfo(ncstatxyt( 8,:),'vpwpxyt'    ,'Turbulent mom. flux'         ,'m^2/s^2','mt'  )
      call ncinfo(ncstatxyt( 9,:),'uwxyt'      ,'Kinematic mom. flux'         ,'m^2/s^2','mt'  )
      call ncinfo(ncstatxyt( 10,:),'wthlxyt'   ,'Kinematic heat flux'         ,'K m/s'  ,'mt'  )
      call ncinfo(ncstatxyt( 11,:),'vwxyt'     ,'Kinematic mom. flux'         ,'m^2/s^2','mt'  )
      call ncinfo(ncstatxyt( 12,:),'usgsxyt'   ,'SGS mom. flux'               ,'m^2/s^2','mt'  )
      call ncinfo(ncstatxyt( 13,:),'thlsgsxyt' ,'SGS heat flux'               ,'K m/s'  ,'mt'  )
      call ncinfo(ncstatxyt( 14,:),'vsgsxyt'   ,'SGS mom. flux'               ,'K m/s'  ,'mt'  )
      call ncinfo(ncstatxyt( 15,:),'thlpthlptxy','Temp. variance'             ,'K^2'    ,'tt'  )
      call ncinfo(ncstatxyt( 16,:),'upuptxyc'  ,'u variance'                  ,'m^2/s^2','tt'  )
      call ncinfo(ncstatxyt( 17,:),'vpvptxyc'  ,'v variance'                  ,'m^2/s^2','tt'  )
      call ncinfo(ncstatxyt( 18,:),'wpwptxyc'  ,'w variance'                  ,'m^2/s^2','tt'  )
      call ncinfo(ncstatxyt( 19,:),'tketxyc'   ,'tke'                         ,'m^2/s^2','tt'  )

      if (myid==0) then      
        call open_nc(xytname, ncidxyt, nrecxyt, n3=khigh-klow+1)
        if (nrecxyt==0) then
          call define_nc( ncidxyt, 1, tncstatxyt)
          call writestat_dims_nc(ncidxyt)
        end if
        call define_nc( ncidxyt, nstatxyt, ncstatxyt)
      end if
    end if

    !> Generate time, y and x averaged NetCDF: zdump.xxx.nc
    if (ltdump) then
 
      tname(7:9) = cmyid
      tname(11:13) = cexpnr
      call ncinfo(tncstatt(1,:),'time'      ,'Time'                        ,'s'      ,'time')
      call ncinfo(ncstatt( 1,:),'ut'        ,'Streamwise velocity'         ,'m/s'    ,'mttt'  )
      call ncinfo(ncstatt( 2,:),'vt'        ,'Spanwise velocity'           ,'m/s'    ,'tmtt'  )
      call ncinfo(ncstatt( 3,:),'wt'        ,'Vertical velocity'           ,'m/s'    ,'ttmt'  )
      call ncinfo(ncstatt( 4,:),'sca1t'     ,'Concentration field 1'       ,'g/m^3'  ,'tttt'  )
      call ncinfo(ncstatt( 5,:),'sca2t'     ,'Concentration field 2'       ,'g/m^3'  ,'tttt'  )
      call ncinfo(ncstatt( 6,:),'sca3t'     ,'Concentration field 3'       ,'g/m^3'  ,'tttt'  )
      call ncinfo(ncstatt( 7,:),'sca4t'     ,'Concentration field 4'       ,'g/m^3'  ,'tttt'  )
      call ncinfo(ncstatt( 8,:),'wpsca4pt'  ,'Turbulent flux 4'            ,'gm/s'   ,'ttmt'  )
      call ncinfo(ncstatt( 9,:),'sv4sgs'    ,'SGS flux 4'                  ,'gm/s'   ,'ttmt'  )
      call ncinfo(ncstatt( 10,:),'wpsca1pt' ,'Turbulent flux 1'            ,'gm/s'   ,'ttmt'  )
      call ncinfo(ncstatt( 11,:),'wpsca2pt' ,'Turbulent flux 2'            ,'gm/s'   ,'ttmt'  )
      call ncinfo(ncstatt( 12,:),'wpsca3pt' ,'Turbulent flux 3'            ,'gm/s'   ,'ttmt'  )
      call ncinfo(ncstatt( 13,:),'sv1sgs'   ,'SGS flux 1'                  ,'gm/s'   ,'ttmt'  )
      call ncinfo(ncstatt( 14,:),'sv2sgs'   ,'SGS flux 2'                  ,'gm/s'   ,'ttmt'  )
      call ncinfo(ncstatt( 15,:),'sv3sgs'   ,'SGS flux 3'                  ,'gm/s'   ,'ttmt'  )

!      if (myid==0) then      
        call open_nc(tname, ncidt, nrect, n1=imax, n2=jmax, n3=khigh-klow+1)
        if (nrect==0) then
          call define_nc( ncidt, 1, tncstatt)
          call writestat_dims_nc(ncidt)
        end if
        call define_nc( ncidt, nstatt, ncstatt)
!      end if
    end if

    !> Generate time, y and x averaged NetCDF for tke budget: tkedump.xxx.nc
    if (ltkedump) then
    
      tkename(9:11) = cexpnr
      call ncinfo(tncstattke(1,:),'time' ,'Time'                                 ,'s'       ,'time')
      call ncinfo(ncstattke( 1,:),'p_b'  ,'p_bant production or consumption term', 'm^2/s^3','tt'  )
      call ncinfo(ncstattke( 2,:),'t_p'  ,'total viscous transport (?)'          , 'm^2/s^3','tt'  )
      call ncinfo(ncstattke( 3,:),'adv'  ,'Advection by mean wind'               , 'm^2/s^3','tt'  )
      call ncinfo(ncstattke( 4,:),'t_t'  ,'Total turb???'                        , 'm^2/s^3','tt'  )
      call ncinfo(ncstattke( 5,:),'t_sgs','total SGS  term'                      , 'm^2/s^3','tt'  )
      call ncinfo(ncstattke( 6,:),'p_t'  ,'Shear production term'                , 'm^2/s^3','tt'  )
      call ncinfo(ncstattke( 7,:),'t_v'  ,'Resolved viscous dissipation term'    , 'm^2/s^3','tt'  )
      call ncinfo(ncstattke( 8,:),'d_sgs','SGS dissipation term'                 , 'm^2/s^3','tt'  )

      if (myid==0) then      
        call open_nc(tkename, ncidtke, nrectke, n3=khigh-klow+1)
        if (nrectke==0) then
          call define_nc( ncidtke, 1, tncstattke)
          call writestat_dims_nc(ncidtke)
        end if
        call define_nc( ncidtke, nstattke, ncstattke)
      endif !myid==0
   
    endif

    if (lslicedump) then

      slicename(11:13) = cmyid
      slicename(15:17) = cexpnr

      call ncinfo(tncstatslice(1,:),'time'     ,'Time'   ,'s'   ,'time')
      call ncinfo(ncstatslice( 1,:),'sca_kb1'  ,'Scalar field at kb', '-', 'tt0t')
      call ncinfo(ncstatslice( 2,:),'sca_ave1' ,'Averaged scalar field over canyon', '-', 'tt0t')
      call ncinfo(ncstatslice( 3,:),'sca_kb2'  ,'Scalar field at kb+1', '-', 'tt0t')
      call ncinfo(ncstatslice( 4,:),'sca_ave2' ,'Averaged scalar field over canyon', '-', 'tt0t')
      call ncinfo(ncstatslice( 5,:),'sca_kb3'  ,'Scalar field at kb+1', '-', 'tt0t')
      call ncinfo(ncstatslice( 6,:),'sca_ave3' ,'Averaged scalar field over canyon', '-', 'tt0t')
      call ncinfo(ncstatslice( 7,:),'u_kb'     ,'Streamwise velocity at kb', '-', 'mt0t')
      call ncinfo(ncstatslice( 8,:),'v_kb'     ,'Spanwise velocity at kb', '-', 'tm0t')

      call open_nc(slicename, ncidslice, nrecslice, n1=imax, n2=jmax)

      if (nrecslice==0) then
        call define_nc( ncidslice, 1, tncstatslice)
        call writestat_dims_nc(ncidslice)  
      end if

      call define_nc( ncidslice, nstatslice, ncstatslice)

    end if

    !> Set times to zero so works for warm starts... could have issues with warmstarts here...
    tsamplep = 0.
    tstatsdumpp = 0.

  end subroutine initstatsdump

  !-------------------------
  !> Generate and write statistics into NetCDF file format
  !-------------------------
 
  subroutine statsdump

  use modfields,        only : um,up,vm,wm,svm,qtm,thlm,pres0,ncstaty,ncstatxy,ncstatyt,ncstattke,&
                               ncstatslice,t_t,t_v,t_p,t_sgs,d_sgs,p_b,p_t,adv,IIc,IIu,IIv,&
                               IIw,IIuw,IIvw,IIct,IIwt,IIut,IIvt,IIuwt,IIcs,IIws,IIus,IIvs,IIuws,&
                               IIvws,slice,slice2,slice3,slice4,slice5,slice6,slice7,slice8,&
                               uyt,vyt,wyt,thlyt,qtyt,&
                               sca1yt,sca2yt,sca3yt,thlsgsyt,usgsyt,&
                               usgsxyt,thlsgsxyt,vsgsxyt,uwtik,vwtjk,wtjk,vtjk,&
                               wthltk,thlthlt,utik,wtik,wmt,thltk,thlt,uxyt,vxyt,wxyt,thlxyt,&
                               ncstatxyt,qtxyt,ncstatt,uutc,vvtc,wwtc,utc,vtc,wtc,&
                               umt,vmt,sv1t,sv2t,sv3t,sv4t,sv1tk,sv2tk,sv3tk,sv4tk,wsv1tk,wsv2tk,wsv3tk,wsv4tk,&
                               sv1sgst,sv2sgst,sv3sgst,sv4sgst
  use modglobal,        only : ib,ie,ih,ihc,xf,xh,jb,je,jhc,jgb,jge,dy,dyi,jh,ke,kb,kh,khc,rk3step,&
                               timee,cexpnr,tsample,tstatsdump,jtot,imax,jmax,dzf,&
                               ltempeq,zh,dxf,dzf,lmassflowr,dzh2i,lprofforc,lscasrcl,&
                               lslicedump,lchem,dzhi,dzhiq,dxhi,massflowrate,lmoist,nsv
!  use modsubgriddata,   only : ekm,sbshr
  use modstat_nc,       only : writestat_nc,writestat_1D_nc
  use modmpi,           only : myid,cmyid,my_real,mpi_sum,avey_ibm,mpierr,&
                               comm3d,avexy_ibm,nprocs
  use modsurfdata,      only : thls
  use modsubgrid,       only : ekh,ekm
  use modstatistics,    only : genstats,tkestats
  implicit none

  !> Create fields to be used in statistics

  ! interpolated fields
!  real, dimension(ib:ie,jb:je,kb:ke)     :: umc
!  real, dimension(ib:ie,jb:je,kb:ke)     :: vmc
!  real, dimension(ib:ie,jb:je,kb:ke)     :: wmc
  real, dimension(ib:ie,jb:je,kb:ke+kh)     :: thlk
  real, dimension(ib:ie,jb:je,kb:ke+kh)     :: uik
  real, dimension(ib:ie,jb:je,kb:ke+kh)     :: wik
  real, dimension(ib:ie,jb:je,kb:ke+kh)     :: vjk
  real, dimension(ib:ie,jb:je,kb:ke+kh)     :: wjk
  real, dimension(ib:ie,jb:je,kb:ke+kh)     :: uc
  real, dimension(ib:ie,jb:je,kb:ke+kh)     :: vc
  real, dimension(ib:ie,jb:je,kb:ke+kh)     :: wc

  ! SGS fluxes
  real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)     :: thlsgs
  real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)     :: usgs
  real, dimension(ib-ih:ie+ih,jb-jh:je+jh,kb:ke+kh)     :: vsgs

  ! t-averaged fields
  real, dimension(ib:ie,jb:je,kb:ke+kh)        :: sv1k
  real, dimension(ib:ie,jb:je,kb:ke+kh)        :: sv2k
  real, dimension(ib:ie,jb:je,kb:ke+kh)        :: sv3k
  real, dimension(ib:ie,jb:je,kb:ke+kh)        :: sv4k
  real, dimension(ib:ie,jb:je,kb:ke+kh)        :: wpsv1p
  real, dimension(ib:ie,jb:je,kb:ke+kh)        :: wpsv2p
  real, dimension(ib:ie,jb:je,kb:ke+kh)        :: wpsv3p
  real, dimension(ib:ie,jb:je,kb:ke+kh)        :: wpsv4p
  real, dimension(ib:ie,jb:je,kb:ke+kh)        :: sv1sgs
  real, dimension(ib:ie,jb:je,kb:ke+kh)        :: sv2sgs
  real, dimension(ib:ie,jb:je,kb:ke+kh)        :: sv3sgs
  real, dimension(ib:ie,jb:je,kb:ke+kh)        :: sv4sgs

  ! y-averaged fields
  real, dimension(ib:ie,kb:ke)                 :: uy
  real, dimension(ib:ie,kb:ke)                 :: vy
  real, dimension(ib:ie,kb:ke)                 :: wy
  real, dimension(ib:ie,kb:ke)                 :: thly
  real, dimension(ib:ie,kb:ke)                 :: qty
  real, dimension(ib:ie,kb:ke)                 :: sca1y
  real, dimension(ib:ie,kb:ke)                 :: sca2y
  real, dimension(ib:ie,kb:ke)                 :: sca3y
  real, dimension(ib:ie,kb:ke)                 :: usgsy
  real, dimension(ib:ie,kb:ke)                 :: thlsgsy

  real, dimension(ib:ie,kb:ke)                 :: uwyik
  real, dimension(ib:ie,kb:ke)                 :: wthlyk
  real, dimension(ib:ie,kb:ke)                 :: wyik
  real, dimension(ib:ie,kb:ke)                 :: uyik
  real, dimension(ib:ie,kb:ke)                 :: thlyk
  real, dimension(ib:ie,kb:ke)                 :: upwpyik
  real, dimension(ib:ie,kb:ke)                 :: wpthlpyk 

  ! xy-averaged fields
  real, dimension(kb:ke+kh)                    :: uxy
  real, dimension(kb:ke+kh)                    :: vxy
  real, dimension(kb:ke+kh)                    :: wxy
  real, dimension(kb:ke+kh)                    :: thlxy
  real, dimension(kb:ke+kh)                    :: qtxy
  real, dimension(kb:ke+kh)                    :: usgsxy
  real, dimension(kb:ke+kh)                    :: thlsgsxy
  real, dimension(kb:ke+kh)                    :: vsgsxy
  real, dimension(kb:ke+kh)                    :: sca1xy

  real, dimension(kb:ke+kh)                    :: uwxyik
  real, dimension(kb:ke+kh)                    :: vwxyjk
  real, dimension(kb:ke+kh)                    :: wthlxyk
  real, dimension(kb:ke+kh)                    :: thlxyk
  real, dimension(kb:ke+kh)                    :: wxyik
  real, dimension(kb:ke+kh)                    :: uxyik
  real, dimension(kb:ke+kh)                    :: vxyjk
  real, dimension(kb:ke+kh)                    :: wxyjk
  real, dimension(kb:ke+kh)                    :: upwpxyik
  real, dimension(kb:ke+kh)                    :: wpthlpxyk
  real, dimension(kb:ke+kh)                    :: vpwpxyjk

  ! fluxes
  real, dimension(ib:ie,kb:ke)                 :: upwptyik
  real, dimension(ib:ie,kb:ke)                 :: wpthlptyk
  real, dimension(ib:ie,kb:ke)                 :: uwtyik
  real, dimension(ib:ie,kb:ke)                 :: wthltyk
  real, dimension(kb:ke+kh)                    :: upwptxyik
  real, dimension(kb:ke+kh)                    :: wpthlptxyk
  real, dimension(kb:ke+kh)                    :: thlpthlptxy
  real, dimension(kb:ke+kh)                    :: upuptxyc
  real, dimension(kb:ke+kh)                    :: vpvptxyc
  real, dimension(kb:ke+kh)                    :: wpwptxyc
  real, dimension(kb:ke+kh)                    :: tketxyc
  real, dimension(kb:ke+kh)                    :: vpwptxyjk
  real, dimension(kb:ke+kh)                    :: uwtxyik
  real, dimension(kb:ke+kh)                    :: wthltxyk
  real, dimension(kb:ke+kh)                    :: vwtxyjk

  real, allocatable :: field(:,:), varsy(:,:,:),varsyt(:,:,:),varstke(:,:),varsxy(:,:),&
                       varslice(:,:,:),varsxyt(:,:),varst(:,:,:,:)
  real    :: tstatsdumppi,emom
  integer :: i,j,k,ip,im,jp,jm,kp,km
  integer :: writecounter = 1
  integer :: reclength

  if (.not. rk3step==3)  return

  if (tsamplep > tsample) then

    if (lytdump .or. lydump .or. lxydump .or. lxytdump) then

      tstatsdumppi = 1./tstatsdumpp 

      !> Perform required interpolations for flux calculations
      !  tg3315 for variable x and z-grids this needs to change
      do k=kb,ke+kh
        do j=jb,je
          do i=ib,ie

            uik(i,j,k) = 0.5*dzhi(k)*(um(i,j,k)*dzf(k-1) + um(i,j,k-1)*dzf(k))
            wik(i,j,k) = 0.5*dxhi(i)*(wm(i,j,k)*dxf(i-1) + wm(i-1,j,k)*dxf(i))
            vjk(i,j,k) = 0.5*dzhi(k)*(vm(i,j,k)*dzf(k-1) + vm(i,j,k-1)*dzf(k))
            wjk(i,j,k) = 0.5*        (wm(i,j,k)          + wm(i,j-1,k))
            uc (i,j,k) = 0.5*dxhi(i)*(um(i,j,k)*dxf(i-1) + um(i-1,j,k)*dxf(i))
            vc (i,j,k) = 0.5*        (vm(i,j,k)          + vm(i,j-1,k))
            wc (i,j,k) = 0.5*dzhi(k)*(wm(i,j,k)*dzf(k-1) + wm(i,j,k-1)*dzf(k))

            ! SGS fluxes
            ! interps ekm to cell corner (uw)
            emom = ( dzf(k-1) * ( ekm(i,j,k)*dxf(i-1)  + ekm(i-1,j,k)*dxf(i) )  + &
                     dzf(k)   * ( ekm(i,j,k-1)*dxf(i-1) + ekm(i-1,j,k-1)*dxf(i) ) )*dxhi(i) * dzhiq(k)
            usgs(i,j,k)  = emom * ( (um(i,j,k)-um(i,j,k-1)) *dzhi(k) &
                        +(wm(i,j,k)-wm(i-1,j,k))  *dxhi(i))

            ! interps ekm to cell corner (vw)
            emom = ( dzf(k-1) * ( ekm(i,j,k)  + ekm(i,j-1,k) )  + &
                     dzf(k)   * ( ekm(i,j,k-1) + ekm(i,j-1,k-1) ) ) * dzhiq(k)

            vsgs(i,j,k)  = emom * ( (vm(i,j,k)-vm(i,j,k-1)) *dzhi(k) &
                        +(wm(i,j,k)-wm(i,j-1,k))  *dyi)

         end do
       end do
     end do

    if (ltempeq) then
      do k=kb,ke+kh
        do j=jb,je
          do i=ib,ie
              thlk(i,j,k) = 0.5*dzhi(k)*(thlm(i,j,k)*dzf(k-1) + thlm(i,j,k-1)*dzf(k))
          end do
        end do
      end do
      do k=kb,ke
        !> SGS fluxes
        thlsgs(:,:,k) = 0.5 * (dzf(k-1)*ekh(:,:,k) + dzf(k)*ekh(:,:,k-1)) &
                        * (thlm(:,:,k)-thlm(:,:,k-1)) * dzh2i(k)
      end do
    end if

    if (nsv>0) then
      do k=kb,ke+kh
        do j=jb,je
          do i=ib,ie
              sv1k(i,j,k) = 0.5*dzhi(k)*(svm(i,j,k,1)*dzf(k-1) + svm(i,j,k-1,1)*dzf(k))
          end do
        end do
      end do
      do k=kb,ke
        sv1sgs(:,:,k) = 0.5 * (dzf(k-1)*ekh(:,:,k) + dzf(k)*ekh(:,:,k-1)) &
                        * (svm(:,:,k,1)-svm(:,:,k-1,1)) * dzh2i(k)
      end do
    end if

    if (nsv>1) then
      do k=kb,ke+kh
        do j=jb,je
          do i=ib,ie
              sv2k(i,j,k) = 0.5*dzhi(k)*(svm(i,j,k,2)*dzf(k-1) + svm(i,j,k-1,2)*dzf(k))
          end do
        end do
      end do
      do k=kb,ke
        sv2sgs(:,:,k) = 0.5 * (dzf(k-1)*ekh(:,:,k) + dzf(k)*ekh(:,:,k-1)) &
                        * (svm(:,:,k,2)-svm(:,:,k-1,2)) * dzh2i(k)
      end do
    end if

    if (nsv>2) then
      do k=kb,ke+kh
        do j=jb,je
          do i=ib,ie
              sv3k(i,j,k) = 0.5*dzhi(k)*(svm(i,j,k,3)*dzf(k-1) + svm(i,j,k-1,3)*dzf(k))
          end do
        end do
      end do
      do k=kb,ke
        sv3sgs(:,:,k) = 0.5 * (dzf(k-1)*ekh(:,:,k) + dzf(k)*ekh(:,:,k-1)) &
                        * (svm(:,:,k,3)-svm(:,:,k-1,3)) * dzh2i(k)
      end do
    end if

    if (nsv>3) then
      do k=kb,ke+kh
        do j=jb,je
          do i=ib,ie
              sv4k(i,j,k) = 0.5*dzhi(k)*(svm(i,j,k,4)*dzf(k-1) + svm(i,j,k-1,4)*dzf(k))
          end do
        end do
      end do
      do k=kb,ke
        sv4sgs(:,:,k) = 0.5 * (dzf(k-1)*ekh(:,:,k) + dzf(k)*ekh(:,:,k-1)) &
                        * (svm(:,:,k,4)-svm(:,:,k-1,4)) * dzh2i(k)
      end do
    end if

      !!>> CALCS FOR INST. STATS
      !> Note: More computationally efficient to spatially average mean quantities first &
      !        for time dependant stats, hence the .or.s.

      !> Average in y-direction      
      if (lydump .or. lytdump) then
     
        uy=0.;vy=0.;wy=0.;thly=0.;qty=0.;sca1y=0.;sca2y=0.;sca3y=0.;thlsgsy=0.;usgsy=0.
 
        call avey_ibm(uy,um(ib:ie,jb:je,kb:ke),ib,ie,jb,je,kb,ke,IIu(ib:ie,jb:je,kb:ke),IIut(ib:ie,kb:ke))
        call avey_ibm(vy,vm(ib:ie,jb:je,kb:ke),ib,ie,jb,je,kb,ke,IIv(ib:ie,jb:je,kb:ke),IIvt(ib:ie,kb:ke))
        call avey_ibm(wy,wm(ib:ie,jb:je,kb:ke),ib,ie,jb,je,kb,ke,IIw(ib:ie,jb:je,kb:ke),IIwt(ib:ie,kb:ke))
        call avey_ibm(thly,thlm(ib:ie,jb:je,kb:ke),ib,ie,jb,je,kb,ke,IIc(ib:ie,jb:je,kb:ke),IIct(ib:ie,kb:ke)) 
        if (lmoist) then
          call avey_ibm(qty,qtm(ib:ie,jb:je,kb:ke),ib,ie,jb,je,kb,ke,IIc(ib:ie,jb:je,kb:ke),IIct(ib:ie,kb:ke))
        end if
        call avey_ibm(uwyik,uik(ib:ie,jb:je,kb:ke)*wik(ib:ie,jb:je,kb:ke),ib,ie,jb,je,kb,ke,IIuw(ib:ie,jb:je,kb:ke),IIuwt(ib:ie,kb:ke))
        call avey_ibm(wthlyk,wm(ib:ie,jb:je,kb:ke)*thlk(ib:ie,jb:je,kb:ke),ib,ie,jb,je,kb,ke,IIw(ib:ie,jb:je,kb:ke),IIwt(ib:ie,kb:ke))
        if (nsv>0) then
          call avey_ibm(sca1y,svm(ib:ie,jb:je,kb:ke,1),ib,ie,jb,je,kb,ke,IIc(ib:ie,jb:je,kb:ke),IIct(ib:ie,kb:ke))
        end if
        if (nsv>1) then
          call avey_ibm(sca2y,svm(ib:ie,jb:je,kb:ke,2),ib,ie,jb,je,kb,ke,IIc(ib:ie,jb:je,kb:ke),IIct(ib:ie,kb:ke))
        end if 
        if (nsv>2) then
          call avey_ibm(sca3y,svm(ib:ie,jb:je,kb:ke,3),ib,ie,jb,je,kb,ke,IIc(ib:ie,jb:je,kb:ke),IIct(ib:ie,kb:ke))
        end if 
        call avey_ibm(usgsy,usgs(ib:ie,jb:je,kb:ke),ib,ie,jb,je,kb,ke,IIuw(ib:ie,jb:je,kb:ke),IIuwt(ib:ie,kb:ke))
        call avey_ibm(thlsgsy,thlsgs(ib:ie,jb:je,kb:ke),ib,ie,jb,je,kb,ke,IIw(ib:ie,jb:je,kb:ke),IIwt(ib:ie,kb:ke))
      end if ! lydump .or. lytdump


      if (lydump) then

        call avey_ibm(uwyik,uik(ib:ie,jb:je,kb:ke)*wik(ib:ie,jb:je,kb:ke),ib,ie,jb,je,kb,ke,IIuw(ib:ie,jb:je,kb:ke),IIuwt(ib:ie,kb:ke))
        call avey_ibm(wthlyk,wm(ib:ie,jb:je,kb:ke)*thlk(ib:ie,jb:je,kb:ke),ib,ie,jb,je,kb,ke,IIw(ib:ie,jb:je,kb:ke),IIwt(ib:ie,kb:ke))
        call avey_ibm(uyik,uik(ib:ie,jb:je,kb:ke),ib,ie,jb,je,kb,ke,IIuw(ib:ie,jb:je,kb:ke),IIuwt(ib:ie,kb:ke))
        call avey_ibm(wyik,wik(ib:ie,jb:je,kb:ke),ib,ie,jb,je,kb,ke,IIuw(ib:ie,jb:je,kb:ke),IIuwt(ib:ie,kb:ke))
        call avey_ibm(thlyk,thlk(ib:ie,jb:je,kb:ke),ib,ie,jb,je,kb,ke,IIw(ib:ie,jb:je,kb:ke),IIwt(ib:ie,kb:ke))

        upwpyik = uwyik - uyik*wyik
        wpthlpyk = wthlyk - wy*thlyk

        where (IIwt==0)
          wpthlpyk  = -999.0                                                                       
        endwhere

        where (IIuwt==0)
          upwpyik    = -999.0
        endwhere

      end if ! lydump

      !> tg3315 10.07.18 - in any case where averaging spatially can assume homogeneity and therefore average
      !  spatially first? Perhaps not due to UCL...? Would save space but goes against triple decomposition 
      !  definition
      !> Average in x and y-direction
      if (lxydump .or. lxytdump) then

        uxy=0.;vxy=0.;wxy=0.;thlxy=0.;qtxy=0.;sca1xy=0.;thlsgsxy=0.;usgsxy=0.;vsgsxy=0.

        !> Spatial averages of mean quantities
        call avexy_ibm(uxy(kb:ke+kh),um(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIu(ib:ie,jb:je,kb:ke+kh),IIus(kb:ke+kh),.true.)
        call avexy_ibm(vxy(kb:ke+kh),vm(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIv(ib:ie,jb:je,kb:ke+kh),IIvs(kb:ke+kh),.true.)
        call avexy_ibm(wxy(kb:ke+kh),wm(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.true.)
        if (ltempeq) then
          call avexy_ibm(thlxy(kb:ke+kh),thlm(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.true.)
        end if
        if (lmoist) then
          call avexy_ibm(qtxy(kb:ke+kh),qtm(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.true.)
        end if
        call avexy_ibm(usgsxy(kb:ke+kh),usgs(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIuw(ib:ie,jb:je,kb:ke+kh),IIuws(kb:ke+kh),.true.)
        call avexy_ibm(thlsgsxy(kb:ke+kh),thlsgs(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.true.)
        call avexy_ibm(vsgsxy(kb:ke+kh),vsgs(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIvw(ib:ie,jb:je,kb:ke+kh),IIvws(kb:ke+kh),.true.)

      end if ! lxydump .or. ztdump

      if (lxydump) then
      
        call avexy_ibm(uwxyik(kb:ke+kh),uik(ib:ie,jb:je,kb:ke+kh)*wik(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIuw(ib:ie,jb:je,kb:ke+kh),IIuws(kb:ke+kh),.true.)
        call avexy_ibm(vwxyjk(kb:ke+kh),vjk(ib:ie,jb:je,kb:ke+kh)*wjk(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIvw(ib:ie,jb:je,kb:ke+kh),IIvws(kb:ke+kh),.true.)
        call avexy_ibm(uxyik(kb:ke+kh),uik(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIuw(ib:ie,jb:je,kb:ke+kh),IIuws(kb:ke+kh),.true.)
        call avexy_ibm(wxyik(kb:ke+kh),wik(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIuw(ib:ie,jb:je,kb:ke+kh),IIuws(kb:ke+kh),.true.)
        call avexy_ibm(wxyjk(kb:ke+kh),wjk(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIvw(ib:ie,jb:je,kb:ke+kh),IIvws(kb:ke+kh),.true.)
        call avexy_ibm(vxyjk(kb:ke+kh),vjk(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIvw(ib:ie,jb:je,kb:ke+kh),IIvws(kb:ke+kh),.true.)
        call avexy_ibm(wthlxyk(kb:ke+kh),wm(ib:ie,jb:je,kb:ke+kh)*thlk(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.true.)
        call avexy_ibm(thlxyk(kb:ke+kh),thlk(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.true.)

        upwpxyik = uwxyik - uxyik*wxyik
        vpwpxyjk = vwxyjk - vxyjk*wxyjk
        wpthlpxyk = wthlxyk - wxy*thlxyk

        where (IIuws==0)
          upwpxyik = -999.
          vpwpxyjk = -999.
        end where

        where (IIws==0)
          wpthlpxyk = -999.
        end where

      end if ! lxydump

      !!>> CALCS FOR TIME DEPENDANT STATS

      !> Average 1-D fields in time
      if (lxytdump) then

        uxyt(kb:ke+kh) = (uxyt(kb:ke+kh)*(tstatsdumpp-tsamplep) + uxy(kb:ke+kh)*tsamplep)*tstatsdumppi
        vxyt(kb:ke+kh) = (vxyt(kb:ke+kh)*(tstatsdumpp-tsamplep) + vxy(kb:ke+kh)*tsamplep)*tstatsdumppi
        wxyt(kb:ke+kh) = (wxyt(kb:ke+kh)*(tstatsdumpp-tsamplep) + wxy(kb:ke+kh)*tsamplep)*tstatsdumppi
        thlxyt(kb:ke+kh) = (thlxyt(kb:ke+kh)*(tstatsdumpp-tsamplep) + thlxy(kb:ke+kh)*tsamplep)*tstatsdumppi
        qtxyt(kb:ke+kh) = (qtxyt(kb:ke+kh)*(tstatsdumpp-tsamplep) + qtxy(kb:ke+kh)*tsamplep)*tstatsdumppi
        usgsxyt(kb:ke+kh) = (usgsxyt(kb:ke+kh)*(tstatsdumpp-tsamplep) + usgsxy(kb:ke+kh)*tsamplep)*tstatsdumppi                                                                                           
        thlsgsxyt(kb:ke+kh) = (thlsgsxyt(kb:ke+kh)*(tstatsdumpp-tsamplep) + thlsgsxy(kb:ke+kh)*tsamplep)*tstatsdumppi
        vsgsxyt(kb:ke+kh) = (vsgsxyt(kb:ke+kh)*(tstatsdumpp-tsamplep) + vsgsxy(kb:ke+kh)*tsamplep)*tstatsdumppi
      end if ! lxytdump

      !> Average 2-D fields in time
      if (lytdump) then

        uyt(ib:ie,kb:ke) = (uyt(ib:ie,kb:ke)*(tstatsdumpp-tsamplep) + uy(ib:ie,kb:ke)*tsamplep)*tstatsdumppi
        vyt(ib:ie,kb:ke) = (vyt(ib:ie,kb:ke)*(tstatsdumpp-tsamplep) + vy(ib:ie,kb:ke)*tsamplep)*tstatsdumppi
        wyt(ib:ie,kb:ke) = (wyt(ib:ie,kb:ke)*(tstatsdumpp-tsamplep) + wy(ib:ie,kb:ke)*tsamplep)*tstatsdumppi
        thlyt(ib:ie,kb:ke) = (thlyt(ib:ie,kb:ke)*(tstatsdumpp-tsamplep) + thly(ib:ie,kb:ke)*tsamplep)*tstatsdumppi
        qtyt(ib:ie,kb:ke) = (qtyt(ib:ie,kb:ke)*(tstatsdumpp-tsamplep) + qty(ib:ie,kb:ke)*tsamplep)*tstatsdumppi
        sca1yt(ib:ie,kb:ke) = (sca1yt(ib:ie,kb:ke)*(tstatsdumpp-tsamplep) + sca1y(ib:ie,kb:ke)*tsamplep)*tstatsdumppi
        sca2yt(ib:ie,kb:ke) = (sca2yt(ib:ie,kb:ke)*(tstatsdumpp-tsamplep) + sca2y(ib:ie,kb:ke)*tsamplep)*tstatsdumppi
        sca3yt(ib:ie,kb:ke) = (sca3yt(ib:ie,kb:ke)*(tstatsdumpp-tsamplep) + sca3y(ib:ie,kb:ke)*tsamplep)*tstatsdumppi
        usgsyt(ib:ie,kb:ke) = (usgsyt(ib:ie,kb:ke)*(tstatsdumpp-tsamplep) + usgsy(ib:ie,kb:ke)*tsamplep)*tstatsdumppi
        thlsgsyt(ib:ie,kb:ke) = (thlsgsyt(ib:ie,kb:ke)*(tstatsdumpp-tsamplep) + thlsgsy(ib:ie,kb:ke)*tsamplep)*tstatsdumppi
      end if !lytdump

      ! Average 3-D fields in time
      if (lxytdump .or. lytdump) then
        uwtik(:,:,kb:ke+kh) = (uwtik(:,:,kb:ke+kh)*(tstatsdumpp-tsamplep) + wik(:,:,kb:ke+kh)*uik(:,:,kb:ke+kh)*tsamplep)*tstatsdumppi
        vwtjk(:,:,kb:ke+kh) = (vwtjk(:,:,kb:ke+kh)*(tstatsdumpp-tsamplep) + wjk(:,:,kb:ke+kh)*vjk(:,:,kb:ke+kh)*tsamplep)*tstatsdumppi
        wthltk(ib:ie,jb:je,kb:ke+kh) = (wthltk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + thlk(ib:ie,jb:je,kb:ke+kh)*wm(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        thlthlt(ib:ie,jb:je,kb:ke+kh) = (thlthlt(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + thlm(ib:ie,jb:je,kb:ke+kh)*thlm(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        uutc(ib:ie,jb:je,kb:ke+kh) = (uutc(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + uc(ib:ie,jb:je,kb:ke+kh)*uc(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        vvtc(ib:ie,jb:je,kb:ke+kh) = (vvtc(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + vc(ib:ie,jb:je,kb:ke+kh)*vc(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        wwtc(ib:ie,jb:je,kb:ke+kh) = (wwtc(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + wc(ib:ie,jb:je,kb:ke+kh)*wc(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        utik(:,:,kb:ke+kh) = (utik(:,:,kb:ke+kh)*(tstatsdumpp-tsamplep) + uik(:,:,kb:ke+kh)*tsamplep)*tstatsdumppi
        wtik(:,:,kb:ke+kh) = (wtik(:,:,kb:ke+kh)*(tstatsdumpp-tsamplep) + wik(:,:,kb:ke+kh)*tsamplep)*tstatsdumppi
        vtjk(:,:,kb:ke+kh) = (vtjk(:,:,kb:ke+kh)*(tstatsdumpp-tsamplep) + vjk(:,:,kb:ke+kh)*tsamplep)*tstatsdumppi
        wtjk(:,:,kb:ke+kh) = (wtjk(:,:,kb:ke+kh)*(tstatsdumpp-tsamplep) + wjk(:,:,kb:ke+kh)*tsamplep)*tstatsdumppi
        wmt(ib:ie,jb:je,kb:ke+kh) = (wmt(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + wm(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        thltk(ib:ie,jb:je,kb:ke+kh) = (thltk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + thlk(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        thlt(ib:ie,jb:je,kb:ke+kh) = (thlt(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + thlm(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        utc(ib:ie,jb:je,kb:ke+kh) = (utc(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + uc(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        vtc(ib:ie,jb:je,kb:ke+kh) = (vtc(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + vc(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        wtc(ib:ie,jb:je,kb:ke+kh) = (wtc(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + wc(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi

      end if !lxytdump .or. lytdump

      ! Other 3-D fields specifically for tdump
      if (ltdump) then
        umt(ib:ie,jb:je,kb:ke+kh) = (umt(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + um(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        vmt(ib:ie,jb:je,kb:ke+kh) = (vmt(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + vm(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        wmt(ib:ie,jb:je,kb:ke+kh) = (wmt(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + wm(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        sv1t(ib:ie,jb:je,kb:ke+kh) = (sv1t(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + svm(ib:ie,jb:je,kb:ke+kh,1)*tsamplep)*tstatsdumppi
        sv2t(ib:ie,jb:je,kb:ke+kh) = (sv2t(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + svm(ib:ie,jb:je,kb:ke+kh,2)*tsamplep)*tstatsdumppi
        sv3t(ib:ie,jb:je,kb:ke+kh) = (sv3t(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + svm(ib:ie,jb:je,kb:ke+kh,3)*tsamplep)*tstatsdumppi
        sv4t(ib:ie,jb:je,kb:ke+kh) = (sv4t(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + svm(ib:ie,jb:je,kb:ke+kh,4)*tsamplep)*tstatsdumppi
        sv1tk(ib:ie,jb:je,kb:ke+kh) = (sv1tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv1k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        sv2tk(ib:ie,jb:je,kb:ke+kh) = (sv2tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv2k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        sv3tk(ib:ie,jb:je,kb:ke+kh) = (sv3tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv3k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        sv4tk(ib:ie,jb:je,kb:ke+kh) = (sv4tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv4k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        wsv1tk(ib:ie,jb:je,kb:ke+kh) = (wsv1tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + wm(ib:ie,jb:je,kb:ke+kh)*sv1k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        wsv2tk(ib:ie,jb:je,kb:ke+kh) = (wsv2tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + wm(ib:ie,jb:je,kb:ke+kh)*sv2k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        wsv3tk(ib:ie,jb:je,kb:ke+kh) = (wsv3tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + wm(ib:ie,jb:je,kb:ke+kh)*sv3k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        wsv4tk(ib:ie,jb:je,kb:ke+kh) = (wsv4tk(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + wm(ib:ie,jb:je,kb:ke+kh)*sv4k(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        sv1sgst(ib:ie,jb:je,kb:ke+kh) = (sv1sgst(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv1sgs(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        sv2sgst(ib:ie,jb:je,kb:ke+kh) = (sv2sgst(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv2sgs(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        sv3sgst(ib:ie,jb:je,kb:ke+kh) = (sv3sgst(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv3sgs(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
        sv4sgst(ib:ie,jb:je,kb:ke+kh) = (sv4sgst(ib:ie,jb:je,kb:ke+kh)*(tstatsdumpp-tsamplep) + sv4sgs(ib:ie,jb:je,kb:ke+kh)*tsamplep)*tstatsdumppi
      end if ! ltdump

      !where (IIuwt==0)
      !  upwpyik    = -999
      !  upwpytik   = -999
      !endwhere
      
      ! EXAMPLE FOR OTHER SLICE PLANES
      !> slice over purifier
!      if (nprocs>7) then
!      if (myid==7) then
!        sca1y(ib:ie,kb:ke) = (sca1y(ib:ie,kb:ke)*(tstatsdumpp-tsamplep) + svm(ib:ie,2,kb:ke,1)*tsamplep)*tstatsdumppi  
!      end if
!      end if

    end if ! lytdump .or. lydump .or. lxydump .or. lxytdump

    ! slicedump fields are generalised so can define what is required here
    if (lslicedump) then
      slice = (slice*(tstatsdumpp-tsamplep) + (0.5*(svm(ib:ie,jb:je,kb,1)+svm(ib:ie,jb:je,kb+1,1)))*tsamplep)*tstatsdumppi
      slice2 = (slice2*(tstatsdumpp-tsamplep) + (sum(svm(ib:ie,jb:je,kb:kb+8,1),3)/9.)*tsamplep)*tstatsdumppi
      slice3 = (slice3*(tstatsdumpp-tsamplep) + (0.5*(svm(ib:ie,jb:je,kb,2)+svm(ib:ie,jb:je,kb+1,2)))*tsamplep)*tstatsdumppi
      slice4 = (slice4*(tstatsdumpp-tsamplep) + (sum(svm(ib:ie,jb:je,kb:kb+8,2),3)/9.)*tsamplep)*tstatsdumppi
      slice5 = (slice5*(tstatsdumpp-tsamplep) + (0.5*(svm(ib:ie,jb:je,kb,3)+svm(ib:ie,jb:je,kb+1,3)))*tsamplep)*tstatsdumppi
      slice6 = (slice6*(tstatsdumpp-tsamplep) + (sum(svm(ib:ie,jb:je,kb:kb+8,3),3)/9.)*tsamplep)*tstatsdumppi
      slice7 = (slice7*(tstatsdumpp-tsamplep) + (um(ib:ie,jb:je,kb)+um(ib:ie,jb:je,kb+1))*tsamplep)*tstatsdumppi
      slice8 = (slice8*(tstatsdumpp-tsamplep) + (vm(ib:ie,jb:je,kb)+vm(ib:ie,jb:je,kb+1))*tsamplep)*tstatsdumppi
    endif !lslicedump

    ! Write y-averaged statistics every tsample
    if (lydump) then
      if (myid == 0) then

        allocate(field(ib:ie,kb:ke)) 
        allocate(varsy(imax,khigh-klow+1,nstaty))

        varsy(:,:,1) = uy(ib:ie,kb:ke)
        varsy(:,:,2) = vy(ib:ie,kb:ke)
        varsy(:,:,3) = wy(ib:ie,kb:ke)
        varsy(:,:,4) = thly(ib:ie,kb:ke)
        varsy(:,:,5) = qty(ib:ie,kb:ke)
        varsy(:,:,6) = sca1y(ib:ie,kb:ke)
        varsy(:,:,7) = sca2y(ib:ie,kb:ke)
        varsy(:,:,8) = sca3y(ib:ie,kb:ke)
        varsy(:,:,9) = upwpyik(ib:ie,kb:ke)
        varsy(:,:,10) = wpthlpyk(ib:ie,kb:ke)
        varsy(:,:,11) = usgsy(ib:ie,kb:ke)
        varsy(:,:,12) = thlsgsy(ib:ie,kb:ke)
        varsy(:,:,13) = uwyik(ib:ie,kb:ke)
        varsy(:,:,14) = wthlyk(ib:ie,kb:ke)

        call writestat_nc(ncidy,1,tncstaty,(/timee/),nrecy,.true.)
        call writestat_nc(ncidy,nstaty,ncstaty,varsy,nrecy,imax,khigh-klow+1)

        deallocate(field,varsy)

      endif !myid
    endif !lydump

    ! Write xy-averaged statistics every tsample
    if (lxydump) then
      if (myid == 0) then
        call writestat_nc(ncidxy,1,tncstatxy,(/timee/),nrecxy,.true.)

        allocate(varsxy(khigh-klow+1,nstatxy))
          varsxy(:,1)  = uxy(kb:ke)
          varsxy(:,2)  = vxy(kb:ke)
          varsxy(:,3)  = wxy(kb:ke)
          varsxy(:,4)  = thlxy(kb:ke)
          varsxy(:,5)  = qtxy(kb:ke)
          varsxy(:,6)  = upwpxyik(kb:ke)
          varsxy(:,7)  = wpthlpxyk(kb:ke)
          varsxy(:,8)  = vpwpxyjk(kb:ke)
          varsxy(:,9)  = usgsxy(kb:ke)
          varsxy(:,10) = thlsgsxy(kb:ke) !wdthldtc(kb:ke)
          varsxy(:,11) = vsgsxy(kb:ke)
          varsxy(:,12) = uwxyik(kb:ke)
          varsxy(:,13) = wthlxyk(kb:ke)
          varsxy(:,14) = vwxyjk(kb:ke)

          call writestat_1D_nc(ncidxy,nstatxy,ncstatxy,varsxy,nrecxy,khigh-klow+1)
      end if !myid
    end if !lxydump

    if (ltkedump) then
      !call genstats(tsamplep,tstatsdumpp,umc,vmc,wmc)
    endif
    tsamplep = dt
  else !timestatsdumpp < tsample

    tsamplep = tsamplep + dt

  endif

  if (tstatsdumpp > tstatsdump) then
    
    ! Final calculations and write xyt-averaged statistics every tsample
    if (lxytdump) then   

      !> Advective flux
      call avexy_ibm(wthltxyk(kb:ke+kh),wmt(ib:ie,jb:je,kb:ke+kh)*thltk(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.true.)
      call avexy_ibm(uwtxyik(kb:ke+kh),utik(ib:ie,jb:je,kb:ke+kh)*wtik(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIuw(ib:ie,jb:je,kb:ke+kh),IIuws(kb:ke+kh),.true.)
      call avexy_ibm(vwtxyjk(kb:ke+kh),vtjk(ib:ie,jb:je,kb:ke+kh)*wtjk(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIvw(ib:ie,jb:je,kb:ke+kh),IIvws(kb:ke+kh),.true.)

      !> Turbulent fluxes
      call avexy_ibm(wpthlptxyk(kb:ke+kh),wthltk(ib:ie,jb:je,kb:ke+kh)-wmt(ib:ie,jb:je,kb:ke+kh)*thltk(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIw(ib:ie,jb:je,kb:ke+kh),IIws(kb:ke+kh),.true.)
      call avexy_ibm(upwptxyik(kb:ke+kh),uwtik(ib:ie,jb:je,kb:ke+kh)-utik(ib:ie,jb:je,kb:ke+kh)*wtik(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIuw(ib:ie,jb:je,kb:ke+kh),IIuws(kb:ke+kh),.true.)
      call avexy_ibm(vpwptxyjk(kb:ke+kh),vwtjk(ib:ie,jb:je,kb:ke+kh)-vtjk(ib:ie,jb:je,kb:ke+kh)*wtjk(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIvw(ib:ie,jb:je,kb:ke+kh),IIvws(kb:ke+kh),.true.)
      call avexy_ibm(thlpthlptxy(kb:ke+kh),thlthlt(ib:ie,jb:je,kb:ke+kh)-thlt(ib:ie,jb:je,kb:ke+kh)*thlt(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.true.)
      call avexy_ibm(upuptxyc(kb:ke+kh),uutc(ib:ie,jb:je,kb:ke+kh)-utc(ib:ie,jb:je,kb:ke+kh)*utc(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.true.)
      call avexy_ibm(vpvptxyc(kb:ke+kh),vvtc(ib:ie,jb:je,kb:ke+kh)-vtc(ib:ie,jb:je,kb:ke+kh)*vtc(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.true.)
      call avexy_ibm(wpwptxyc(kb:ke+kh),wwtc(ib:ie,jb:je,kb:ke+kh)-wtc(ib:ie,jb:je,kb:ke+kh)*wtc(ib:ie,jb:je,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.true.)
      call avexy_ibm(tketxyc(kb:ke+kh),0.5*((wwtc(ib:ie,jb:je,kb:ke+kh)-wtc(ib:ie,jb:je,kb:ke+kh)*wtc(ib:ie,jb:je,kb:ke+kh))+(vvtc(ib:ie,jb:je,kb:ke+kh)-vtc(ib:ie,jb:je,kb:ke+kh)*vtc(ib:ie,jb:je,kb:ke+kh))+(uutc(ib:ie,jb:je,kb:ke+kh)-utc(ib:ie,jb:je,kb:ke+kh)*utc(ib:ie,jb:je,kb:ke+kh))),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc(ib:ie,jb:je,kb:ke+kh),IIcs(kb:ke+kh),.true.)


      if (myid == 0) then
        call writestat_nc(ncidxyt,1,tncstatxyt,(/timee/),nrecxyt,.true.)

        allocate(varsxyt(khigh-klow+1,nstatxyt))
          varsxyt(:,1)  = uxyt(kb:ke)
          varsxyt(:,2)  = vxyt(kb:ke)
          varsxyt(:,3)  = wxyt(kb:ke)
          varsxyt(:,4)  = thlxyt(kb:ke)
          varsxyt(:,5)  = qtxyt(kb:ke)
          varsxyt(:,6)  = upwptxyik(kb:ke)
          varsxyt(:,7)  = wpthlptxyk(kb:ke)
          varsxyt(:,8)  = vpwptxyjk(kb:ke)
          varsxyt(:,9)  = uwtxyik(kb:ke)
          varsxyt(:,10) = wthltxyk(kb:ke) !wdthldtc(kb:ke)
          varsxyt(:,11) = vwtxyjk(kb:ke)
          varsxyt(:,12) = usgsxyt(kb:ke) !wdthldtw(kb:ke)
          varsxyt(:,13) = thlsgsxyt(kb:ke)
          varsxyt(:,14) = vsgsxyt(kb:ke)
          varsxyt(:,15) = thlpthlptxy(kb:ke)
          varsxyt(:,16) = upuptxyc(kb:ke)
          varsxyt(:,17) = vpvptxyc(kb:ke)
          varsxyt(:,18) = wpwptxyc(kb:ke)
          varsxyt(:,19) = tketxyc(kb:ke)
          call writestat_1D_nc(ncidxyt,nstatxyt,ncstatxyt,varsxyt,nrecxyt,khigh-klow+1)
      end if !myid
    end if !lxytdump

    ! Final calculations and write yt-averaged statistics every tsample
    if (lytdump) then

!    call MPI_BCAST(sca1yt ,(ke+kh-(kb-kh))*(ie+ih-(ib-ih)),MY_REAL   ,7,comm3d,mpierr)

      ! Turbulent flux
      call avey_ibm(upwptyik,uwtik(ib:ie,jb:je,kb:ke)-utik(ib:ie,jb:je,kb:ke)*wtik(ib:ie,jb:je,kb:ke),ib,ie,jb,je,kb,ke,IIuw(ib:ie,jb:je,kb:ke),IIuwt(ib:ie,kb:ke))
      call avey_ibm(wpthlptyk,wthltk(ib:ie,jb:je,kb:ke)-wmt(ib:ie,jb:je,kb:ke)*thltk(ib:ie,jb:je,kb:ke),ib,ie,jb,je,kb,ke,IIw(ib:ie,jb:je,kb:ke),IIwt(ib:ie,kb:ke))

      ! Advective flux
      call avey_ibm(uwtyik,utik(ib:ie,jb:je,kb:ke)*wtik(ib:ie,jb:je,kb:ke),ib,ie,jb,je,kb,ke,IIuw(ib:ie,jb:je,kb:ke),IIuwt(ib:ie,kb:ke))
      call avey_ibm(wthltyk,wmt(ib:ie,jb:je,kb:ke)*thltk(ib:ie,jb:je,kb:ke),ib,ie,jb,je,kb,ke,IIw(ib:ie,jb:je,kb:ke),IIwt(ib:ie,kb:ke))

      if (myid == 0) then
          allocate(varsyt(imax,khigh-klow+1,nstaty))
          call writestat_nc(ncidyt,1,tncstatyt,(/timee/),nrecyt,.true.)
          varsyt(:,:,1)  = uyt(ib:ie,kb:ke)
          varsyt(:,:,2)  = vyt(ib:ie,kb:ke)
          varsyt(:,:,3)  = wyt(ib:ie,kb:ke)
          varsyt(:,:,4)  = thlyt(ib:ie,kb:ke)
          varsyt(:,:,5)  = qtyt(ib:ie,kb:ke)
          varsyt(:,:,6)  = sca1yt(ib:ie,kb:ke)
          varsyt(:,:,7)  = sca2yt(ib:ie,kb:ke)
          varsyt(:,:,8)  = sca3yt(ib:ie,kb:ke)
          varsyt(:,:,9)  = upwptyik(ib:ie,kb:ke)
          varsyt(:,:,10) = wpthlptyk(ib:ie,kb:ke)
          varsyt(:,:,11) = uwtyik(ib:ie,kb:ke)
          varsyt(:,:,12) = wthltyk(ib:ie,kb:ke)
          varsyt(:,:,13) = usgsyt(ib:ie,kb:ke)
          varsyt(:,:,14) = thlsgsyt(ib:ie,kb:ke)
          call writestat_nc(ncidyt,nstatyt,ncstatyt,varsyt,nrecyt,imax,khigh-klow+1)
        end if !myid
      end if !lytdump

    ! Final calculations and write t-averaged statistics every tsample
    if (ltdump) then

    wpsv1p = wsv1tk - wmt*sv1tk
    wpsv2p = wsv2tk - wmt*sv2tk
    wpsv3p = wsv3tk - wmt*sv3tk
    wpsv4p = wsv4tk - wmt*sv4tk

!      if (myid == 0) then
          allocate(varst(imax,jmax,khigh-klow+1,nstatt))
          call writestat_nc(ncidt,1,tncstatt,(/timee/),nrect,.true.)
          varst(:,:,:,1)  = umt(ib:ie,jb:je,kb:ke)
          varst(:,:,:,2)  = vmt(ib:ie,jb:je,kb:ke)
          varst(:,:,:,3)  = wmt(ib:ie,jb:je,kb:ke)
          varst(:,:,:,4)  = sv1t(ib:ie,jb:je,kb:ke)
          varst(:,:,:,5)  = sv2t(ib:ie,jb:je,kb:ke)
          varst(:,:,:,6)  = sv3t(ib:ie,jb:je,kb:ke)
          varst(:,:,:,7)  = sv4t(ib:ie,jb:je,kb:ke)
          varst(:,:,:,8)  = wpsv4p(ib:ie,jb:je,kb:ke)
          varst(:,:,:,9)  = sv4sgst(ib:ie,jb:je,kb:ke)
          varst(:,:,:,10) = wpsv1p(ib:ie,jb:je,kb:ke)
          varst(:,:,:,11) = wpsv2p(ib:ie,jb:je,kb:ke)
          varst(:,:,:,12) = wpsv3p(ib:ie,jb:je,kb:ke)
          varst(:,:,:,13) = sv1sgst(ib:ie,jb:je,kb:ke)
          varst(:,:,:,14) = sv2sgst(ib:ie,jb:je,kb:ke)
          varst(:,:,:,15) = sv3sgst(ib:ie,jb:je,kb:ke)
          call writestat_nc(ncidt,nstatt,ncstatt,varst,nrect,imax,jmax,khigh-klow+1)
!        end if !myid
         deallocate(varst)
      end if !ltdump

      if (ltkedump) then
        call tkestatsdump
        if (myid == 0) then
          call writestat_nc(ncidtke,1,tncstattke,(/timee/),nrectke,.true.)
          allocate(varstke(khigh-klow+1,nstattke))
          varstke(:,1) = p_b(kb:ke+kh)
          varstke(:,2) = t_p(kb:ke+kh)
          varstke(:,3) = adv(kb:ke+kh)
          varstke(:,4) = t_t(kb:ke+kh)
          varstke(:,5) = t_sgs(kb:ke+kh)
          varstke(:,6) = p_t(kb:ke+kh)
          varstke(:,7) = t_v(kb:ke+kh)
          varstke(:,8) = d_sgs(kb:ke+kh)              
          call writestat_1D_nc(ncidtke,nstattke,ncstattke,varstke,nrectke,khigh-klow+1)
        end if !myid
      endif !ltkedump 

      if (lslicedump) then

        allocate(varslice(imax,jmax,nstatslice))
        call writestat_nc(ncidslice,1,tncstatslice,(/timee/),nrecslice,.true.)

        varslice(:,:,1) = slice
        varslice(:,:,2) = slice2
        varslice(:,:,3) = slice3
        varslice(:,:,4) = slice4
        varslice(:,:,5) = slice5   
        varslice(:,:,6) = slice6   
        varslice(:,:,7) = slice7   
        varslice(:,:,8) = slice8

!        write(*,*), myid
!        write(*,*), 'ncidslice,1,tncstatslice,(/timee/),nrecslice,.true.', ncidslice,1,tncstatslice,(/timee/),nrecslice

        call writestat_nc(ncidslice,nstatslice,ncstatslice,varslice,nrecslice,imax,jmax)

!        deallocate(varslice)

      endif

      tstatsdumpp = dt

    else !tstatsdumpp < tstatsdump
   
      tstatsdumpp = tstatsdumpp + dt
  
    endif  

  end subroutine statsdump
   
  !> tg3315 still under going work to be completed
  subroutine tkestatsdump

  use modfields,        only : u0,v0,w0,thl0,uav,vav,wav,uuav,vvav,wwav,uvav,uwav,vwav,thlav,thlthlav,pres0,thluav,thlvav,thlwav,&
                               upupav,vpvpav,wpwpav,thlpthlpav,upvpav,upwpav,vpwpav,thlpupav,thlpvpav,thlpwpav,presav,&
                               strain2av,disssgsav,t_vav,tvmx,tvmy,tvmz,tsgsmx1,tsgsmx2,tsgsmy1,tsgsmy2,tsgsmz1,t_sgsav,nusgsav,&
                               tpm,t_pav,ttmx,ttmy,ttmz,t_tav,p_bav,d_sgsav,p_tav,tkeadv,tsgsmz1,tsgsmz2,t_t,t_v,t_p,t_sgs,d_sgs,&
                               p_b,p_t,adv,IIc,IIcs
  use modglobal,        only : ib,ie,ih,jb,je,jgb,jge,dy,jh,ke,kb,kh,rk3step,timee,cexpnr,tsample,tstatsdump,jtot,imax,dzf,&
                               dzf,dzfi,dzhi,dxf,dxfi,dyi,dxhi,dy2i,grav,numol
  use modmpi,           only : myid,cmyid,my_real,mpi_sum,avey_ibm,mpierr,comm3d,excjs,avexy_ibm
  use modsurfdata,      only : thls
  use modsubgrid,       only : ekh
  implicit none

  real, dimension(ib:ie,jb:je,kb:ke)  :: disssgsfl     ! average subgrid visc. * average rate of strain squared : 2*<nu_t>*<Sij>*<Sij>
  real, dimension(ib:ie,jb:je,kb:ke)  :: dissresav     ! average resolved dissipation: 2*nu*<Sij'*Sij'> = 2*nu*( <Sij*Sij> - <Sij>*<Sij> )
    real, dimension(ib:ie,jb:je,kb:ke)  :: tke           ! tke = 0.5*<ui'ui'>
    real, dimension(ib:ie,jb:je,kb:ke)  :: mke           ! = <ui>d/dxj(<ui><uj>) + <ui>d/dxj(<ui'uj'>) = <ui>d/dxj(<ui*uj>)
    real, dimension(ib:ie+1,jb  :je,  kb:ke)    :: dummyx
    real, dimension(ib:ie,  jb-1:je+1,kb:ke)    :: dummyy
    real, dimension(ib:ie,  jb  :je,  kb:ke+1)  :: dummyz

  integer i,j,k,ip,im,jp,jm,kp,km
  real strainav2
  real dummy

    ! Tvav = (Tvm - <ui>*d/dxj(<Sij>)  ) + 2*nu*<Sij'Sij'>
    ! Tvm = Tvmx + Tvmy + Tvmz -> therefore: subtraction, then interpolation,
    ! then addition of 2*nu*<Sij'Sij'>
    do k=kb,ke
      km = k-1
      kp = k+1
      do j=jb,je
        jp = j+1
        jm = j-1
        do i=ib,ie
          im = i-1
          ip = i+1

!            t_vav(i,j,k) =  0.5*( (tvmx(i,j,k) - (                      &
             dummyx(i,j,k) =  (                      &
                              ( numol  * (uav(i+1,j,k)-uav(i,j,k))*dxfi(i) &
                              -numol * (uav(i,j,k)-uav(i-1,j,k))*dxfi(i-1) ) * 2. * dxhi(i) &
                              + &
                              ( numol * ( (uav(i,jp,k)-uav(i,j,k))   *dyi &
                              +(vav(i,jp,k)-vav(i-1,jp,k))*dxhi(i)) &
                              - numol * ( (uav(i,j,k)-uav(i,jm,k))   *dyi &
                              +(vav(i,j,k)-vav(i-1,j,k))  *dxhi(i)) &
                              ) * dyi &
                              + &
                              ( numol * ( (uav(i,j,kp)-uav(i,j,k))   *dzhi(kp) &
                              +(wav(i,j,kp)-wav(i-1,j,kp))*dxhi(i)) &
                              - numol * ( (uav(i,j,k)-uav(i,j,km))   *dzhi(k) &
                              +(wav(i,j,k)-wav(i-1,j,k))  *dxhi(i)) &
                              ) *dzfi(k) )
               ! y-direction
               dummyy(i,j,k) =  (                        &
                                ( numol * ( (vav(i+1,j,k)-vav(i,j,k))   *dxhi(i+1) &
                                +(uav(i+1,j,k)-uav(i+1,jm,k))*dyi) &
                                -numol * ( (vav(i,j,k)-vav(i-1,j,k))   *dxhi(i) &
                                +(uav(i,j,k)-uav(i,jm,k))    *dyi) &
                                ) * dxfi(i) &        ! = d/dx( Km*(dv/dx + du/dy) )
                                + &
                                (numol * (vav(i,jp,k)-vav(i,j,k)) &
                                -numol * (vav(i,j,k)-vav(i,jm,k))  ) * 2. * dy2i &        ! =d/dy( 2*Km*(dv/dy) )
                                + &
                                ( numol * ( (vav(i,j,kp)-vav(i,j,k))    *dzhi(kp) &
                                +(wav(i,j,kp)-wav(i,jm,kp))  *dyi) &
                                -numol * ( (vav(i,j,k)-vav(i,j,km))    *dzhi(k) &
                                +(wav(i,j,k)-wav(i,jm,k))    *dyi)   &
                                ) * dzfi(k) )                    ! = d/dz( Km*(dv/dz + dw/dy) )
               ! z-direction
               dummyz(i,j,k) = (                        &
                              ( numol * ( (wav(i+1,j,k)-wav(i,j,k))    *dxhi(i+1) &
                              +(uav(i+1,j,k)-uav(i+1,j,km)) *dzhi(k) ) &
                              -numol * ( (wav(i,j,k)-wav(i-1,j,k))    *dxhi(i) &
                              +(uav(i,j,k)-uav(i,j,km))     *dzhi(k) ) &
                             )*dxfi(i) &
                             + &
                             ( numol * ( (wav(i,jp,k)-wav(i,j,k))     *dyi &
                             +(vav(i,jp,k)-vav(i,jp,km))   *dzhi(k) ) &
                             -numol * ( (wav(i,j,k)-wav(i,jm,k))     *dyi &
                             +(vav(i,j,k)-vav(i,j,km))     *dzhi(k) ) &
                             )*dyi &
                             + &
                             ( numol * (wav(i,j,kp)-wav(i,j,k)) *dzfi(k) &
                             -numol * (wav(i,j,k)-wav(i,j,km)) *dzfi(km) ) * 2. &
                             * dzhi(k))

               strainav2 =  ( &
                            ((uav(ip,j,k)-uav(i,j,k))    *dxfi(i)     )**2    + &
                            ((vav(i,jp,k)-vav(i,j,k))    *dyi         )**2    + &
                            ((wav(i,j,kp)-wav(i,j,k))    *dzfi(k)     )**2    )

               strainav2 = strainav2 + 0.125 * ( &
                            ((wav(i,j,kp)-wav(im,j,kp))   *dxhi(i)     + &
                            (uav(i,j,kp)-uav(i,j,k))      *dzhi(kp)  )**2    + &
                            ((wav(i,j,k)-wav(im,j,k))     *dxhi(i)     + &
                            (uav(i,j,k)-uav(i,j,km))      *dzhi(k)   )**2    + &
                            ((wav(ip,j,k)-wav(i,j,k))     *dxhi(ip)     + &
                            (uav(ip,j,k)-uav(ip,j,km))    *dzhi(k)   )**2    + &
                            ((wav(ip,j,kp)-wav(i,j,kp))   *dxhi(ip)     + &
                            (uav(ip,j,kp)-uav(ip,j,k))    *dzhi(kp)  )**2    )

               strainav2 = strainav2 + 0.125 * ( &
                            ((uav(i,jp,k)-uav(i,j,k))     *dyi     + &
                            (vav(i,jp,k)-vav(im,jp,k))    *dxhi(i)        )**2    + &
                            ((uav(i,j,k)-uav(i,jm,k))     *dyi     + &
                            (vav(i,j,k)-vav(im,j,k))      *dxhi(i)        )**2    + &
                            ((uav(ip,j,k)-uav(ip,jm,k))   *dyi     + &
                            (vav(ip,j,k)-vav(i,j,k))      *dxhi(ip)       )**2    + &
                            ((uav(ip,jp,k)-uav(ip,j,k))   *dyi     + &
                            (vav(ip,jp,k)-vav(i,jp,k))    *dxhi(ip)       )**2    )

               strainav2 = strainav2 + 0.125 * ( &
                           ((vav(i,j,kp)-vav(i,j,k))    *dzhi(kp) + &
                           (wav(i,j,kp)-wav(i,jm,kp))   *dyi        )**2    + &
                           ((vav(i,j,k)-vav(i,j,km))    *dzhi(k)+ &
                           (wav(i,j,k)-wav(i,jm,k))     *dyi        )**2    + &
                           ((vav(i,jp,k)-vav(i,jp,km))  *dzhi(k)+ &
                           (wav(i,jp,k)-wav(i,j,k))     *dyi        )**2    + &
                           ((vav(i,jp,kp)-vav(i,jp,k))  *dzhi(kp) + &
                           (wav(i,jp,kp)-wav(i,j,kp))   *dyi        )**2    )

               dissresav(i,j,k) = 2.*numol  *(strain2av(i,j,k) - strainav2)  !resolved dissipation

          end do
        end do
      end do

      ! BC's 
      tvmx   (ie+1,:,:) =  tvmx   (ie,:,:)
      tsgsmx1(ie+1,:,:) =  tsgsmx1(ie,:,:)
      tsgsmx2(ie+1,:,:) =  tsgsmx2(ie,:,:)
      dummyx (ie+1,:,:) =  dummyx (ie,:,:)
      ttmx   (ie+1,:,:) =  ttmx   (ie,:,:)
      call excjs( tvmy   , ib,ie,jb,je,kb,ke,0,1)   ! jb-1 is not used
      call excjs( tsgsmy1, ib,ie,jb,je,kb,ke,0,1)   ! jb-1 is not used
      call excjs( tsgsmy2, ib,ie,jb,je,kb,ke,0,1)   ! jb-1 is not used
      call excjs( dummyy,  ib,ie,jb,je,kb,ke,0,1)   ! jb-1 is not used
      call excjs( ttmy   , ib,ie,jb,je,kb,ke,0,1)   ! jb-1 is not used
      tvmz   (:,:,ke+1) =  tvmz   (:,:,ke)
      tsgsmz1(:,:,ke+1) =  tsgsmz1(:,:,ke)
      tsgsmz2(:,:,ke+1) =  tsgsmz2(:,:,ke)
      dummyz (:,:,ke+1) =  dummyz(:,:,ke)
      ttmz   (:,:,ke+1) =  ttmz   (:,:,ke)

      do k=kb,ke
        km = k-1
        kp = k+1
        do j=jb,je
          jp = j+1
          jm = j-1
            do i=ib,ie
             im = i-1
             ip = i+1

             ! Total viscous dissipation
             t_vav(i,j,k) =  0.5*( (tvmx(i, j,k) - dummyx(i,j,k) *uav(i, j,k))  + &
                             (tvmx(ip,j,k) - dummyx(ip,j,k)*uav(ip,j,k))) &
                             + 0.5*( (tvmy(i,j, k) - dummyy(i,j,k) *vav(i,j, k))  + &
                             (tvmy(i,jp,k) - dummyy(i,jp,k)*vav(i,jp,k))) &
                             + 0.5*( (tvmz(i,j,k ) - dummyz(i,j,k) *wav(i,j,k ))  + &
                             (tvmz(i,j,kp) - dummyz(i,j,kp)*wav(i,j,kp))) &
                             + dissresav(i,j,k)         ! d/dxj(2*nu*<ui'Sij'>) = <u_i*d/dxj(2*nu*Sij')> +2*nu*<Sij'Sij'>

!      Now the same for subgrid stress        
!      <d/dxj(2*u_i'*nu_t*Sij)'> = <u_i'*d/dxj(2*nu_t*Sij)'> + <(2*nu_t*Sij)'*Sij'>
!                                = <u_i*d/dxj(2*nu_t*Sij)> -
!                                  <u_i>*<d/dxj(2*nu_t*Sij)>   + <2*nu_t*Sij*Sij> -
!                                  <(2*nu_t*Sij)>*<Sij>
!                                = <u_i*d/dxj(2*nu_t*Sij)> -
!                                  <u_i>*<d/dxj(2*nu_t*Sij)>   + <2*nu_t*Sij*Sij> -
!                                  2*<nu_t>*<Sij>*<Sij> - 2*<nu_t'*Sij'>*<Sij>


     !---------------------------------------
     !Total subgrid TKE
     !---------------------------------------

             ! Mean SGS dissipation
             disssgsfl(i,j,k) = 2.*nusgsav(i,j,k)*strainav2 ! = 2*<nu_sgs>*<sij>*<sij>

              
             ! TKE
             tke(i,j,k)       = 0.5*(0.5*(upupav(ip,j,k)+upupav(i,j,k)) + &
                                  0.5*(vpvpav(i,jp,k)+vpvpav(i,j,k)) + &
                                  0.5*(wpwpav(i,j,kp)+wpwpav(i,j,k)))

             ! total SGS
             t_sgsav(i,j,k) =  0.5*( (tsgsmx1(i,j,k) -  uav(i,j,k) *tsgsmx2(i,j,k)) + &
                              (tsgsmx1(ip,j,k) - uav(ip,j,k)*tsgsmx2(ip,j,k))) &
                               + & ! = <2*nu_t*SijSij> - <2*nu_t*Sij>*<Sij>
                              0.5*( (tsgsmy1(i,j,k) -  vav(i,j,k) *tsgsmy2(i,j,k)) + &
                              (tsgsmy1(i,jp,k) - vav(i,jp,k)*tsgsmy2(i,jp,k))) &
                              + & ! = <2*nu_t*SijSij> - <2*nu_t*Sij>*<Sij>        
                              0.5*( (tsgsmz1(i,j,k) -  vav(i,j,k) *tsgsmz2(i,j,k)) + &
                              (tsgsmz1(i,j,kp) - vav(i,j,kp)*tsgsmz2(i,j,kp))) &
                              + disssgsav(i,j,k) - disssgsfl(i,j,k)
             ! -2*<nu_t'Sij'>*<Sij>  should still be added!
             
             ! SGS dissipation
             d_sgsav(i,j,k)= - disssgsav(i,j,k) + disssgsfl(i,j,k)
             ! +2*<nu_t'Sij'>*<Sij>  should still be added! (is compensated with above)


     !---------------------------------------
     !Total pressure TKE
     !---------------------------------------

             ! Pressure correlation term
             ! - <uj'*dp'/dxj> = - <uj*dp/dxj> + <uj>*d<p>/dxj
             t_pav(i,j,k)   = tpm(i,j,k) + &
                              0.5*(uav(i,j,k)*(presav(i,j,k)-presav(i-1,j,k))*dxhi(i) + &
                              uav(i+1,j,k)*(presav(i+1,j,k)-presav(i,j,k))*dxhi(i+1)) &
                              + &
                              0.5*(vav(i,j,k)*(presav(i,j,k)-presav(i,j-1,k))*dyi + &
                              vav(i,j+1,k)*(presav(i,j+1,k)-presav(i,j,k))*dyi) &
                              + &
                              0.5*(wav(i,j,k)*(presav(i,j,k)-presav(i,j,k-1))*dzhi(k) + &
                              wav(i,j,k+1)*(presav(i,j,k+1)-presav(i,j,k))*dzhi(k+1))
                              ! - d/dxj(<0.5*ui'ui'uj'>) = -<uj'd/dxj(<0.5*ui'ui'>) + <ui'uj'><Sij>
!                             = -<uj*d/dxj(0.5*ui'ui')> + <uj>*d/dxj(<0.5*ui'ui'> +
!                             <ui'uj'><Sij>) 


!            ttav(i,j,k)   = ttm(i,j,k) - 


     !---------------------------------------
     !Total advection TKE
     !---------------------------------------

!            <advection term N.S. times ui> = MKE + A - Pshear - Tt
!            Tt = -<ui'd/dxj(ui'uj')> = -<d/dxj(0.5*ui'ui'uj')> = A + MKE - Pshear - Total

             !Pshear =Ptav = -<ui'uj'>d/dxj(<Sij>) = -<ui'uj'>d<ui>/dxj

             ! mechanical or shear production
             p_tav(i,j,k)    = - ( &
                               0.5 *(upupav(i,j,k)+upupav(ip,j,k))* (uav(ip,j,k)-uav(i,j,k))*dxfi(i)  + & ! <u'u'>*d<u>/dx     
                               0.25*(upvpav(i,j,k)  *(uav(i, j, k)-uav(i, jm,k) )*dyi + &
                               upvpav(i,jp,k) *(uav(i, jp,k)-uav(i, j, k) )*dyi + &
                               upvpav(ip,j,k) *(uav(ip,j, k)-uav(ip,jm,k) )*dyi + &
                               upvpav(ip,jp,k)*(uav(ip,jp,k)-uav(ip,j, k) )*dyi) + & ! <u'v'>*d<u>/dy
                               0.25*(upwpav(i, j,k ) *(uav(i, j,k )-uav(i,j,km))*dzhi(k) + &
                               upwpav(i, j,kp) *(uav(i, j,kp)-uav(i, j,k))*dzhi(kp) + &
                               upwpav(ip,j,k ) *(uav(ip,j,k )-uav(ip,j,km))*dzhi(k) + &
                               upwpav(ip,j,kp) *(uav(ip,j,kp)-uav(ip,j,k))*dzhi(kp)) + & ! <u'w'>*d<u>/dz
                               0.25*(upvpav(i, j, k) *(vav(i, j, k)-vav(im,j,k))*dxhi(i) + &
                               upvpav(ip,j, k) *(vav(ip,j, k)-vav(i, j,k))*dxhi(ip) + &
                               upvpav(i, jp,k) *(vav(i,jp,k)-vav(im,jp,k))*dxhi(i) + &
                               upvpav(ip,jp,k) *(vav(ip,jp,k)-vav(i,jp,k))*dxhi(ip)) + & ! <u'v'>*d<v>/dx
                               0.5 *(vpvpav(i,j,k)+vpvpav(i,jp,k))*(vav(i,jp,k)-vav(i,j,k))*dyi + & ! <v'v'>*d<v>/dy
                               0.5 *(vpvpav(i,j,k)+vpvpav(i,jp,k))*(vav(i,jp,k)-vav(i,j,k))*dyi + & ! <v'v'>*d<v>/dy
                               0.25*(vpwpav(i,j ,k ) *(vav(i,j ,k )-vav(i,j,km))*dzhi(k) + &
                               vpwpav(i,j ,kp) *(vav(i,j ,kp)-vav(i,j ,k))*dzhi(kp) + &
                               vpwpav(i,jp,k ) *(vav(i,jp,k)-vav(i,jp,km))*dzhi(k) + &
                               vpwpav(i,jp,kp) *(vav(i,jp,kp)-vav(i,jp,k))*dzhi(kp)) + & ! <v'w'>*d<v>/dz
                               0.25*(upwpav(i, j, k) *(wav(i, j,k )-wav(im,j,k))*dxhi(i) + &
                               upwpav(ip,j, k) *(wav(ip,j,k )-wav(i, j,k))*dxhi(ip) + &
                               upwpav(i, j,kp) *(wav(i,j,kp)-wav(im,j,kp))*dxhi(i) + &
                               upwpav(ip,j,kp) *(wav(ip,j,kp)-wav(i,j,kp))*dxhi(ip)) + & ! <u'w'>*d<w>/dx
                               0.25*(vpwpav(i,j,k)  *(wav(i,j, k )-wav(i,jm,k ) )*dyi + &
                               vpwpav(i,jp,k) *(wav(i,jp,k )-wav(i,j, k ) )*dyi + &
                               vpwpav(ip,j,k) *(wav(i,j, kp)-wav(i,jm,kp) )*dyi + &
                               vpwpav(ip,jp,k)*(wav(i,jp,kp)-wav(i,j, kp) )*dyi) + & ! <v'w'>*d<w>/dy
                               0.5 *(wpwpav(i,j,k)+wpwpav(i,j,kp))*(wav(i,j,kp)-wav(i,j,k))*dzfi(k) ) ! <w'w'>*d<w>/dz 

             ! Mean kinetic energy term (expected to be small).
             mke(i,j,k)      = 0.5*(uav(ip,j,k)+uav(i,j,k))*(uuav(ip,j,k)-uuav(i,j,k))*dxfi(i)  +        & !<u>*d<uu>/dx          
                               0.5*(uav(i, j,k)*(uvav(i ,jp,k)-uvav(i ,j,k))*dyi  + & ! <u>*d<uv>/dy                      
                               uav(ip,j,k)*(uvav(ip,jp,k)-uvav(ip,j,k))*dyi) +        &
                               0.5*(uav(i, j,k)*(uwav(i ,j,kp)-uwav(i ,j,k))*dzfi(k)  + & ! <u>*d<uw>/dz     
                               uav(ip,j,k)*(uwav(ip,j,kp)-uwav(ip,j,k))*dzfi(k)) +        &
                               0.5*(vav(i,j, k)*(uvav(ip,j ,k)-uvav(i,j ,k))*dxfi(i) + & ! <v>*d<uv>/dx
                               vav(i,jp,k)*(uvav(ip,jp,k)-uvav(i,jp,k))*dxfi(i)) + &
                               0.5*(vav(i,jp,k)+vav(i,j,k))*(vvav(i,jp,k)-vvav(i,j,k))*dyi +        & ! <v>*d<vv>/dy
                               0.5*(vav(i,j ,k)*(vwav(i,j ,kp)-vwav(i,j ,k))*dzfi(k)  + & ! <v>*d<vw>/dz                      
                               vav(i,jp,k)*(vwav(i,jp,kp)-vwav(i,jp,k))*dzfi(k)) + &
                               0.5*(wav(i,j,k )*(uwav(ip,j,k )-uwav(i,j,k ))*dxfi(i) + & ! <w>*d<uw>/dx
                               wav(i,j,kp)*(uwav(ip,j,kp)-uwav(i,j,kp))*dxfi(i)) +        &
                               0.5*(wav(i,j,k )*(vwav(i,jp,k )-vwav(i,j,k ))*dyi  + & ! <w>*d<vw>/dy                      
                               wav(i,j,kp)*(vwav(i,jp,kp)-vwav(i,j,kp))*dyi) +        &
                               0.5*(wav(i,j,kp)+wav(i,j,k))*(wwav(i,j,kp)-wwav(i,j,k))*dzfi(k) ! <w>*d<ww>/dz  

             ! Advection of TKE
             tkeadv(i,j,k)   = 0.5*(uav(i, j,k)*(tke(i, j,k)-tke(im,j,k))*dxhi(i) + & ! <u>*de/dx
                               uav(ip,j,k)*(tke(ip,j,k)-tke(i ,j,k))*dxhi(ip)) +            & !
                               0.5*(vav(i, j,k)*(tke(i,j ,k)-tke(i,jm,k))*dyi     + & ! <v>*de/dy
                               vav(i,jp,k)*(tke(i,jp,k)-tke(i,j ,k))*dyi) +            &
                               0.5*(wav(i,j,k )*(tke(i,j,k )-tke(i,j,km))*dzhi(k) + & ! <w>*de/dz
                               wav(i,j,kp)*(tke(i,j,kp)-tke(i,j,k ))*dzhi(kp))
        
             ! <advection term N.S. times ui> = MKE + A - Pshear - Tt
             ! Tt = -<ui'd/dxj(ui'uj')> = -<d/dxj(0.5*ui'ui'uj')> = A      +    MKE   -
             ! Pshear  -   Total
             !                                                    = tkeadv +    mke   -
             !                                                    p_tav   -   ttm
             !        t_tav(i,j,k)   = tkeadv(i,j,k) + mke(i,j,k) - p_tav(i,j,k) - ttm(i,j,k)                 
            
             t_tav(i,j,k)   = tkeadv(i,j,k) + mke(i,j,k) - p_tav(i,j,k)  &
                              - 0.5*(ttmx(i,j,k) + ttmx(ip,j,k))        &
                              - 0.5*(ttmy(i,j,k) + ttmy(i,jp,k))        &
                              - 0.5*(ttmz(i,j,k) + ttmz(i,j,kp))

             p_bav(i,j,k)   = (grav/thls)*0.5*(thlpwpav(i,j,k)+thlpwpav(i,j,kp)) !use of thls here...????

          end do
        end do
      end do    

    ! need updating tg3315
    call avexy_ibm(p_b(kb:ke+kh),p_bav(:,:,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc,IIcs,.true.)
    call avexy_ibm(t_p(kb:ke+kh),t_pav(:,:,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc,IIcs,.true.)
    call avexy_ibm(adv(kb:ke+kh),tkeadv(:,:,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc,IIcs,.true.)
    call avexy_ibm(t_t(kb:ke+kh),t_tav(:,:,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc,IIcs,.true.)
    call avexy_ibm(t_sgs(kb:ke+kh),t_sgsav(:,:,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc,IIcs,.true.)
    call avexy_ibm(p_t(kb:ke+kh),p_tav(:,:,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc,IIcs,.true.)
    call avexy_ibm(d_sgs(kb:ke+kh),d_sgsav(:,:,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc,IIcs,.true.)
     call avexy_ibm(t_v(kb:ke+kh),t_vav(:,:,kb:ke+kh),ib,ie,jb,je,kb,ke,ih,jh,kh,IIc,IIcs,.true.)

   end subroutine tkestatsdump

  !-------------------------
  !> Clean up when leaving the run
  !------------------------

  subroutine exitstatsdump
      use modstat_nc, only : exitstat_nc
      use modglobal, only  : lslicedump,ltdump
    implicit none

!       if (lydump) then 
!         call exitstat_nc(ncid)
!       endif
 
! will doing this ruin the averaging? ... try tg3315      
!       if (lytdump) then              
!         call exitstat_nc(ncidt)
!       endif

!      if (ltkedump) then              
!        call exitstat_nc(ncidtke)
!      endif

!      if (lslicedump) then              
!        call exitstat_nc(ncidslice)
!      endif

!       if (ltdump) then              
!         call exitstat_nc(ncidt)
!       endif

  end subroutine exitstatsdump


end module modstatsdump
