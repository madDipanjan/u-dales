!> \file modetrees.f90
!!tg3315, ns4513, 20 Mar 2017 

!> Input trees into DALES model.

module modheatpump
implicit none
save

contains

    subroutine heatpump
    use modglobal,  only : ib,ie,jb,je,kb,ke,ih,jh,xh,xf,yh,yf,tree,&
                           nheatpumps,lheatpumps,itot,jtot,zf,zh,&
                           ltempeq, rk3step
    use modfields,  only : wp,thlp
    use modmpi,     only : myid, myidx,myidy,mpi_sum,mpierr,comm3d,mpierr,my_real,nprocx,nprocy
    
    implicit none
    integer :: i,j,k,n,m,il,iu,jl,ju,kl,ku

    if (lheatpumps .eqv. .false.) return
    
    ! if (rk3step==3) then
    ! ! write(*,*) myid, " ", myidx, " ", myidy, " ", ib, " ", ie, " ", jb, " ", je, " ", kb, " ", ke
    ! ! write(*,*) myid, " ", myidx, " ", myidy, " ", xh(ib), " ", xh(ie), " ", yh(jb), " ", yh(je), " ", zh(kb), " ", zh(ke)
    ! ! write(*,*) myid, " ", myidx, " ", myidy, " ", xf(ib), " ", xf(ie), " ", yf(jb), " ", yf(je), " ", zf(kb), " ", zf(ke)
    !     write(*,*) myid, " ", myidx, " ", myidy, " ", xf(ib+myidx*itot/nprocx), " ", xf(ie+myidx*itot/nprocx), &
    !     " ", yf(jb+myidy*jtot/nprocy), " ", yf(je+myidy*jtot/nprocy), " ", zf(kb), " ", zf(ke)
    ! end if

    if (rk3step==3) then
        do k = kb, ke
        do j = jb, je
            ju = j+myidy*jtot/nprocy
        do i = ib, ie
            iu = i+myidx*itot/nprocx
            ! if (k==8 .and. ju==43 .and. iu==97) then
            !     write(*,*) myid, " ", myidx, " ", myidy, " ", i, " ", j, " ", k, " ", &
            !             xf(iu), " ", yf(ju), " ", zf(k)
            ! elseif (k==8 .and. ju==17 .and. iu==58) then
            !     write(*,*) myid, " ", myidx, " ", myidy, " ", i, " ", j, " ", k, " ", &
            !             xf(iu), " ", yf(ju), " ", zf(k)
            ! elseif (k==8 .and. ju==39 .and. iu==32) then
            !     write(*,*) myid, " ", myidx, " ", myidy, " ", i, " ", j, " ", k, " ", &
            !             xf(iu), " ", yf(ju), " ", zf(k)
            ! elseif (k==8 .and. ju==25 .and. iu==82) then
            !     write(*,*) myid, " ", myidx, " ", myidy, " ", i, " ", j, " ", k, " ", &
            !             xf(iu), " ", yf(ju), " ", zf(k)

            ! if ( (k==6) .and. ( (ju==55 .and. iu==97) .or. (ju==54 .and. iu==98) .or. &
            !                     (ju==55 .and. iu==98) .or. (ju==54 .and. iu==99) .or. &
            !                     (ju==66 .and. iu==93) .or. (ju==66 .and. iu==94) .or. &
            !                     (ju==65 .and. iu==94) .or. (ju==65 .and. iu==95) ) )  then
            if ( (k==6) .and. ( (ju==14 .and. iu==49) .or. (ju==17 .and. iu==47) ) )  then
                
                    ! write(*,*) myid, " ", myidx, " ", myidy, " ", i, " ", j, " ", k, " ", &
                    !            xf(iu), " ", yf(ju), " ", zf(k)

                    wp(iu,ju,k) = 0.25
                    ! thlp(iu,ju,k) = 273
                    thlp(i,j,k) = thlp(i,j,k) - 15
            end if
        end do
        end do
        end do
    end if

    end subroutine heatpump

  end module modheatpump
