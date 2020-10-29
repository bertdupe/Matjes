module m_MCstep
use m_derived_types, only : lattice
use m_basic_types, only : vec_point
use m_sampling
use m_choose_spin
use m_relaxtyp
use m_createspinfile
use m_H_public
    implicit none

private
public :: MCstep
contains
!
! ===============================================================
!
SUBROUTINE MCstep(lat,N_spin,E_total,E,Magnetization,kt,acc,rate,nb,cone,ising,equi,overrel,sphere,underrel,Hams)
use m_constants, only : k_b,pi
use mtprng
Implicit none
! input
type(lattice),intent(inout)    :: lat
logical, intent(in) :: equi,overrel,sphere,underrel,ising
real(kind=8), intent(in) :: kt
integer, intent(in) :: N_spin
real(kind=8), intent(inout) :: E_total,Magnetization(3),E(8),acc,rate,cone,nb
class(t_H), intent(in) :: Hams(:)

! internal variable
type(mtprng_state) :: state
!     Energy difference Delta E and -DE/kT, Dmag and Dq
real(kind=8) :: DE,E_new,E_old,Dmag(3)
real(kind=8) :: tmp
!     memory value for Theta
real(kind=8) :: choice
!     flipped Spin
real(kind=8) :: S_new(3),S_old(3)
integer :: iomp

call choose_spin(iomp,N_spin)

! cone angle update and maximal magnetic moment change
if ((rate.gt.0.5d0).and.(cone.lt.pi(1.0d0)))then
     cone=cone+0.0001d0
elseif ((rate.lt.0.50d0).and.(cone.gt.0.01d0)) then
     cone=cone-0.0001d0
endif
!---------------------------------------
! here are the different sampling
! first the sphere sampling
!---------------------------------------
!PB: the way with different if's seems arther risky for choosing the S_new
S_new   =0.0d0
if(lat%nmag>1) ERROR STOP "WILL NOT WORK FOR nmag>1"
if (sphere) then
   S_new=sphereft(lat%M%modes_v(:,iomp),cone)
!---------------------------------------
! Algorithm square sampling
!---------------------------------------
elseif (equi) then
   S_new=equirep(lat%M%modes_v(:,iomp),cone)
endif
!---------------------------------
! different relaxation process
!---------------------------------
if (underrel) then
   S_new=underrelax(iomp,lat)
elseif (overrel) then
   S_new=overrelax(iomp,lat)
endif
if (ising) S_new=-lat%M%modes_v(:,iomp)

!----------------------------------
!       Calculate the energy difference if this was flipped
!       and decider, if the Spin flip will be performed
!----------------------------------
!Energy of old configuration
E_old=energy_single(Hams,iomp,lat)
S_old=lat%M%modes_v(:,iomp)
if(lat%nmag>1) ERROR STOP "THIS WILL NOT WORK WITH NMAG>1"

!Energy of the new configuration
lat%M%modes_v(:,iomp)=S_new
E_new=energy_single(Hams,iomp,lat)
if(lat%nmag>1) ERROR STOP "THIS WILL NOT WORK WITH NMAG>1"

lat%M%modes_v(:,iomp)=S_old
!! variation of the energy for this step
DE=E_new-E_old


nb=nb+1.0d0
if ( accept(kt,DE) ) then
    Dmag=-lat%M%modes_v(:,iomp)+S_new
    if(lat%nmag>1) ERROR STOP "THIS and next line WILL NOT WORK WITH NMAG>1"
! update the spin
    lat%M%modes_v(:,iomp)=S_new
! update the quantities
    E_total=E_total+DE
    Magnetization=Magnetization+Dmag
    acc=acc+1.0d0
endif

rate=acc/nb

END SUBROUTINE MCstep
! ===============================================================


function accept(kt,DE)
    real(8),intent(in)  :: kt,DE

    real(8) :: choice, tmp
    logical :: accept
    
    accept=.False.
    !accept energy gain
    if(dE<0.0d0)then
        accept=.true.
        return
    endif
    ! security in case kt is 0
    if(kt<1.0d-10)then
       return
    endif
    
    tmp=-DE/kT
    !prevent exp(tmp) underflow
    if(tmp<-200.0d0) return
#ifdef CPP_MRG
    choice=mtprng_rand_real1(state)
#else
    CALL RANDOM_NUMBER(Choice)
#endif
    if (Choice.lt.exp(tmp)) Then
        accept=.True.
    endif

end function


end module
