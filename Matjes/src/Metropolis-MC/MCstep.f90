module m_MCstep
use m_derived_types, only : lattice
use m_basic_types, only : vec_point
use m_sampling
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
SUBROUTINE MCstep(lat,N_spin,E_total,E,Magnetization,kt,acc,rate,nb,cone,ising,equi,overrel,sphere,underrel,Hams)
!#ifdef CPP_MPI
!      use m_parameters, only : d_mu,n_ghost,i_separate,i_average,i_ghost
!      use m_make_box, only : ghost_border,rank_nei
!      use m_mpi_prop, only : isize,MPI_COMM,N,start,MPI_COMM_BOX,irank_working,irank_box
!#endif
use m_constants, only : k_b,pi
use mtprng
Implicit none
! input
type(lattice),intent(inout)    :: lat
!type(vec_point), intent(inout) :: mode(:)
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
! decomposition of the energy
real(kind=8) :: E_dec_new(8),E_dec_old(8)
integer :: iomp
!temporary assignment
type(vec_point), pointer :: mode(:)

!#ifdef CPP_MPI
!real(kind=8) :: mpi_DE(isize),mpi_S_store(3,isize),mpi_mu(isize)
!integer :: imin(1),irank_discuss,irank_send
!! shall the comm occur
!logical :: discuss,at_border,discuss_trans
!
!! control of the point to point comm
!integer :: reqs(4)
!
!include 'mpif.h'
!
!mpi_DE=0.0d0
!mpi_S_store=0.0d0
!mpi_mu=0.0d0
!discuss=.False.
!at_border=.False.
!discuss_trans=.False.
!irank_discuss=0
!irank_send=0
!tmp=0.0d0
!
!#endif

mode=>lat%ordpar%all_l_modes

E_dec_new=0.0d0
E_dec_old=0.0d0

!#ifdef CPP_MPI
!call choose_spin(Ilat,i_separate,i_average,i_ghost,i_stone,n_world, &
!    &   mu_s,irank_working,shape_spin,MPI_COMM,spin,start)
!#else
call choose_spin(iomp,N_spin)
!#endif

!     Initialize Store variables
S_new   =0.0d0

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
    !TODO UPDATE
   S_new=underrelax(iomp,lat)
elseif (overrel) then
    !TODO UPDATE
   S_new=overrelax(iomp,lat)
endif

if (ising) S_new=-lat%M%modes_v(:,iomp)
!----------------------------------
!       Calculate the energy difference if this was flipped
!       and decider, if the Spin flip will be performed
!----------------------------------
!Energy of old configuration
E_old=energy_single(Hams,iomp,lat)

!Energy of the new configuration
S_old=mode(iomp)%w(1:3)

mode(iomp)%w(1:3)=S_new

E_new=energy_single(Hams,iomp,lat)

mode(iomp)%w(1:3)=S_old

!! variation of the energy for this step
DE=E_new-E_old

!---------------------------------------------------
! Here is the different way to compute de next 
! candidate for the spin
!---------------------------------------------------

! End of the broadcast
!#ifdef CPP_MPI
!if ((.not.i_separate).and.(.not.i_average).and.(.not.i_ghost)) then
!
!   call mpi_gather(DE,1,MPI_REAL8,mpi_DE,1,MPI_REAL8,0,MPI_COMM,ierr)
!   call mpi_gather(S_new,3,MPI_REAL8,mpi_S_store(:,:),3,MPI_REAL8,0,MPI_COMM,ierr)
!
!   if (irank_working.eq.0) then
!      imin=minloc(mpi_DE)
!      DE=mpi_DE(imin(1))
!      S_new=mpi_S_store(:,imin(1))
!   endif
!
!   call mpi_bcast(DE,1,MPI_REAL8,0,MPI_COMM,ierr)
!   call mpi_bcast(S_new,3,MPI_REAL8,0,MPI_COMM,ierr)
!
!endif
!#endif

nb=nb+1.0d0

if ( accept(kt,DE) ) then

   Dmag=-mode(iomp)%w(1:3)+S_new

! update the spin
   mode(iomp)%w(1:3)=S_new

! update the quantities
   E_total=E_total+DE
   E=E+E_dec_new-E_dec_old
   Magnetization=Magnetization+Dmag
   acc=acc+1.0d0
endif

!#ifdef CPP_MPI
!! if the move is accepted, one has to make sure that the updated spin does not lie at one interface
!
!if (i_ghost) then
!   at_border=ghost_border(1,Ilat(1),Ilat(2),Ilat(3))
!!         write(*,*) at_border ,Ilat(1),Ilat(2), irank_working, irank_box
!   do i=1,n_ghost
!      if ((accept).and.(at_border).and.((i-1).eq.irank_box)) then
!         discuss_trans=.True.
!         irank_send=irank_box
!      endif
!      call mpi_allreduce(discuss_trans,discuss,1,mpi_logical,MPI_LOR,MPI_COMM_BOX,ierr)
!      call mpi_allreduce(irank_send,irank_discuss,1,mpi_integer,MPI_sum,MPI_COMM_BOX,ierr)
!
!! if it lies at one of the interface then update the spins of the neighbors
!
!      if (discuss) call MC_transfer(Ilat,S_new,irank_discuss,spin,shape_spin)
!
!      discuss_trans=.False.
!      irank_send=0
!   enddo
!endif
!#endif

rate=acc/nb
nullify(mode)

!      write(*,*) acc,nb,kt/k_B,irank_working
! part that does the parallel tempering

END SUBROUTINE MCstep
! ===============================================================



function accept(kt,DE)
implicit none
real(kind=8) :: kt,DE
real(kind=8) :: choice, tmp
logical :: accept

accept=.False.

! security in case kt is 0
if (kt.gt.1.0d-10) then
!       Calculate the probability of this flip
   tmp=-DE/kT
   if (DE.lt.0.0d0) tmp=2.0d0
else
   tmp=-100.0d0
!       Accept the change of the new energy is lower than the old one
   if (DE.lt.0.0d0) tmp=2.0d0
endif

if (exp(tmp).gt.1.0d0) then
   accept=.True.
else
!#ifdef CPP_MPI
!#ifdef CPP_MRG
!   choice=mtprng_rand_real1(state)
!#else
!   CALL RANDOM_NUMBER(Choice)
!#endif
!   if ((.not.i_separate).and.(.not.i_average).and.(.not.i_ghost)) call mpi_bcast(Choice,1,MPI_REAL8,0,MPI_COMM,ierr)
!#else
#ifdef CPP_MRG
   choice=mtprng_rand_real1(state)
#else
   CALL RANDOM_NUMBER(Choice)
#endif
!#endif
   if (Choice.lt.exp(tmp)) Then
      accept=.True.
   endif
endif

end function



end module
