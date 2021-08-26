module m_relaxation
use m_derived_types, only : lattice
use m_mc_track_val, only: track_val
use m_hamiltonian_collection, only: hamiltonian
use m_topo_commons
use m_constants, only : k_b
implicit none

private
public :: Relaxation

type :: t_relax
    integer :: Nstep=0      !number of MCsteps done
    real(8) :: E=0.0d0      !energy per site
    real(8) :: M(3)=0.0d0   !magnetization average
    real(8) :: q(2)=0.0d0   !topological charge (positive,negative contibutions)
    real(8) :: rate=0.0d0   !acceptance parameter
    real(8) :: cone=0.0d0   !acceptance parameter
contains 
    procedure   :: set => set_relax
end type

contains
!
! ===============================================================
SUBROUTINE Relaxation(lat,io_MC,N_spin,state_prop,kt,H,Q_neigh,work)
    !Relaxes the magnetic state for a given temperature by calling the MC_step io_MC%N_relaxation*io_MC%T_relax*N_spin times.
    !Intermediate relaxation information can be obtained with io_MC%print_relax  at io_MC%n_sizerelax states.
    use m_Corre
    use m_io_files_utils
    use m_convert
    use m_MCstep
    use m_MC_io,only: MC_input
    use m_work_ham_single, only: work_ham_single
    ! input
    type(lattice),intent(inout)         :: lat
    real(kind=8), intent(in)            :: kT
    type(track_val),intent(inout)       :: state_prop
    type(MC_input),intent(in)           :: io_MC 
    integer, intent(in)                 :: N_spin       !number of spins considered
    type(hamiltonian),intent(inout)     :: H
    integer,intent(in),allocatable      :: Q_neigh(:,:) 
    type(work_ham_single),intent(inout) :: work
    integer         :: N_MCStep                     !number of MCsteps within each outer relaxation loop
    type(T_relax)   :: relax(io_MC%n_sizerelax)     !type to store state at different relaxation steps
    integer         :: n_w_step                     !number of steps for each relaxation progress write
    integer         :: i_relaxation,i_MC,i_store    !loop variables
    
    write(6,'(/,a,f8.4,2x,a,/)') 'starting relaxation for T= ',kT/k_B,'Kelvin'
    
    n_w_step=io_MC%N_relaxation/io_MC%n_sizerelax
    N_MCStep=max(1,io_MC%T_relax*N_spin)  !In case T_relax set to zero at least one MCstep is done
    i_store=0
    
    do i_relaxation=1,io_MC%N_relaxation
        Do i_MC=1,N_MCStep
            Call MCStep(lat,io_MC,N_spin,state_prop,kt,H,work)
        enddo
    
        ! Save information about the Relaxation for the Equilibrium files
        if (io_MC%print_relax) then
            if (mod(i_relaxation,n_w_step).eq.0)then
                i_store=i_store+1
                if(i_store>size(relax)) ERROR STOP "nonsensical combination of io_MC%N_thousand and io_MC%n_sizerelax"
                Call relax(i_store)%set(lat,state_prop,Q_neigh,i_relaxation*N_MCStep)
            endif
        endif
    enddo   ! enddo over the relaxation loop
    
    write(6,'(/,a,f8.4,2x,a,/)') 'System is relaxed for T= ',kT/k_B,'Kelvin'
    
    if (io_MC%print_relax)  Call print_relax_arr(Relax,kT)  !print the Equilibrium files
end subroutine

subroutine set_relax(this,lat,state_prop,Q_neigh,Nstep)
    class(t_relax),intent(inout)    :: this
    type(lattice),intent(in)        :: lat
    type(track_val),intent(in)      :: state_prop
    integer,intent(in),allocatable  :: Q_neigh(:,:)
    integer,intent(in),optional     :: Nstep

    real(8)     :: dumy(5)
    
    this%E=state_prop%E_total
    this%M=state_prop%magnetization/lat%nmag/lat%Ncell

    dumy=0.0d0
    if(allocated(Q_neigh))then
        dumy=get_charge(lat,Q_neigh)
    endif
    this%q=dumy(1:2)
    this%rate=state_prop%rate
    this%cone=state_prop%cone
    if(present(Nstep)) this%Nstep = Nstep
end subroutine

subroutine print_relax_arr(relax,kT)
    use m_convert
    type(t_relax),intent(in)    :: relax(:)
    real(8),intent(in)          :: kT 

    integer :: i, io

    open(newunit=io,file=convert('Equilibriumphi-',kT/k_B,'.dat'))
    write(io,'(A)') "#   1:#steps    2:energy        3:M_x           4:M_y           5:M_z           6:M             7:topo          8:topo+         9:topo-         10: rate        11:cone" 
    do i=1,size(relax)
        write(io,'(I12)'     ,advance='no') relax(i)%Nstep
        write(io,'(E16.8)'  ,advance='no') relax(i)%E
        write(io,'(3E16.8)' ,advance='no') relax(i)%M
        write(io,'(E16.8)'  ,advance='no') norm2(relax(i)%M)
        write(io,'(E16.8)'  ,advance='no') sum(relax(i)%q)
        write(io,'(2E16.8)' ,advance='no') relax(i)%q
        write(io,'(E16.8)'  ,advance='no') relax(i)%rate
        write(io,'(E16.8)'  ,advance='no') relax(i)%cone
        write(io,'(/)'  ,advance='no')
    enddo
    close(io)

end subroutine
end module
