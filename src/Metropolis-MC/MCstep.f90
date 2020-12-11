module m_MCstep
use m_derived_types, only : lattice
use m_mc_track_val, only: track_val
use m_sampling, only: sphereft,equirep
use m_choose_spin
use m_relaxtyp
use m_createspinfile
use m_H_public

private
public :: MCstep
contains
!
! ===============================================================
!
SUBROUTINE MCstep(lat,io_MC,N_spin,state_prop,kt,Hams)
    use m_constants, only : k_b,pi
    use m_MC_io,only: MC_input
    use mtprng
    Implicit none
    ! input
    type(lattice),intent(inout)     :: lat
    type(MC_input),intent(in)       :: io_MC 
    real(kind=8), intent(in)        :: kt
    integer, intent(in)             :: N_spin
    type(track_val),intent(inout)   :: state_prop
    class(t_H), intent(in)          :: Hams(:)
    
    ! internal variable
    !     Energy difference Delta E and -DE/kT, Dmag and Dq
    real(kind=8) :: DE,E_new,E_old,Dmag(3)
    !     flipped Spin
    real(kind=8) :: S_new(3),S_old(3)
    integer :: i_spin   !chosen spin index (1:N_cell*nmag)
    integer :: i_site   !unit cell index of chosen spin index
    
    call choose_spin(i_spin,N_spin)
    i_site=(i_spin-1)/lat%nmag+1
#ifdef CPP_DEBUG
    if(lat%nmag>1) ERROR STOP "A lot of things will not work in this routine with nmag>1"
#endif
    !---------------------------------------
    ! here are the different sampling
    ! first the sphere sampling
    !---------------------------------------
    if(io_MC%ising)then
        S_new=-lat%M%modes_3(:,i_spin)
    elseif(io_MC%underrelax)then
        S_new=underrelax(i_spin,lat,Hams)
    elseif(io_MC%overrelax)then
        S_new=overrelax(i_spin,lat,Hams)
    elseif(io_MC%equi)then
        Call cone_update(state_prop%cone,state_prop%rate)
        S_new=equirep(lat%M%modes_3(:,i_spin),state_prop%cone)
    elseif(io_MC%sphere)then
        Call cone_update(state_prop%cone,state_prop%rate)
        S_new=sphereft(lat%M%modes_3(:,i_spin),state_prop%cone)
    else
        STOP "Invalid MC sampling method"
    endif
    
    !----------------------------------
    !       Calculate the energy difference if this was flipped
    !       and decider, if the Spin flip will be performed
    !----------------------------------
    !Energy of old configuration
    E_old=energy_single(Hams,i_site,lat)
    S_old=lat%M%modes_3(:,i_spin)
    
    !Energy of the new configuration
    lat%M%modes_3(:,i_spin)=S_new
    E_new=energy_single(Hams,i_site,lat)
    lat%M%modes_3(:,i_spin)=S_old
    !! variation of the energy for this step
    DE=E_new-E_old
    
    state_prop%nb=state_prop%nb+1.0d0
    if ( accept(kt,DE) ) then
        Dmag=-S_old+S_new
    ! update the spin
        lat%M%modes_3(:,i_spin)=S_new
    ! update the quantities
        state_prop%E_total=state_prop%E_total+DE
        state_prop%Magnetization=state_prop%Magnetization+Dmag
        state_prop%acc=state_prop%acc+1.0d0
    endif
   
    state_prop%rate=state_prop%acc/state_prop%nb
END SUBROUTINE MCstep
! ===============================================================

function accept(kt,DE)
    use m_get_random, only: get_rand_classic
    real(8),intent(in)  :: kt,DE

    real(8) :: choice, tmp
    logical :: accept
    
#ifdef CPP_DEBUG
    ! security in case kt is 0
    if(kt<1.0d-10)then
        ERROR STOP "kt is too small"
    endif
#endif

    accept=dE<0.0d0     !accept all energy gains
    if(.not.accept)then
        !check with temperature if energy loss is accepted
        tmp=-DE/kT
        tmp=max(tmp,-200.0d0) !prevent exp(tmp) underflow
        choice=get_rand_classic()
        accept=Choice.lt.exp(tmp)
    endif
end function

pure subroutine  cone_update(cone,rate)
    use m_constants, only : pi
    real(8),intent(inout)   :: cone
    real(8),intent(in)      :: rate

    ! cone angle update and maximal magnetic moment change
    if ((rate.gt.0.5d0).and.(cone.lt.pi))then
         cone=cone+0.0001d0
    elseif ((rate.lt.0.50d0).and.(cone.gt.0.01d0)) then
         cone=cone-0.0001d0
    endif
end subroutine

end module
