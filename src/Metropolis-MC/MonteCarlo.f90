module m_montecarlo
implicit none
contains
!
! Routine that does the Monte Carlo (and not the parallel tempering)
!
!subroutine MonteCarlo(my_lattices,mag_motif,io_simu,gra_topo,ext_param)
subroutine montecarlo(my_lattice,io_simu,ext_param,Hams,com)
    use, intrinsic :: iso_fortran_env, only : output_unit
    use mpi_util 
    use m_constants, only : k_b,pi
    use m_derived_types, only : lattice,t_cell,io_parameter,simulation_parameters
    use m_rw_MC
    use m_topo_commons, only : neighbor_Q,get_charge
    use m_MCstep
    use m_H_public
    use m_relaxation
    use m_write_config
    use m_set_temp
    use m_input_types,only: MC_input
    use m_mc_exp_val
    use m_mc_track_val,only: track_val
    use m_average_MC, only: get_neighbours

    type(lattice), intent(inout)            :: my_lattice
    type(io_parameter), intent(in)          :: io_simu
    type(simulation_parameters), intent(in) :: ext_param
    class(t_H), intent(in)                  :: Hams(:)
    class(mpi_type),intent(in)              :: com
    
    !!!!!!!!!!!!!!!!!!!!!!!
    ! internal variables
    !!!!!!!!!!!!!!!!!!!!!!!
    ! slope of the MC
    integer                     :: i_relax,i_kT,i_MC,N_cell,size_table,i_measure
    ! variable for the temperature
    real(8)                     :: kT,kTini,kTfin
    ! variables that being followed during the simulation
    real(8)                     :: qeulerp,qeulerm,dumy(5)
    type(track_val)             :: state_prop
    real(8),pointer :: M3(:,:)
    ! contribution of the different energy parts
    type(exp_val),allocatable   :: measure(:)
    real(8),allocatable         :: kt_all(:)
    type(MC_input)              :: io_MC
    integer,allocatable         :: Q_neigh(:,:)
    integer                     :: filen_kt_acc(2)
    integer                     :: ierr
    integer                     :: iT_global(2),iT_local(2)
    integer,parameter           :: io_status=output_unit
    !mpi parameters
    integer                     :: cnt(com%Np)
    integer                     :: displ(com%Np)
	integer, allocatable :: flat_nei(:,:)
    integer, allocatable :: indexNN(:)
    
    ! initialize the variables
    call rw_MC(io_MC)
    size_table=io_MC%n_Tsteps
    write(io_status,'(/,a,I6,a,/)') "you are calculating",size_table," temperatures"
    allocate(kt_all(size_table),source=0.0d0)
    kTini=ext_param%ktini
    kTfin=ext_param%ktfin
    N_cell=my_lattice%Ncell
    filen_kt_acc=5
    
    ! allocate and fill flat table of first neighbours
	call get_neighbours(my_lattice,flat_nei,indexNN)

    ! initialize the temperatures
    call ini_temp(kt_all,kTfin,kTini,size_table,io_simu%io_warning)
    filen_kt_acc(1)=max(int(log10(maxval(kt_all))),1)
    
    Call neighbor_Q(my_lattice,Q_neigh)
    
    !intialize tracking variables (total energy, magnetization sum)
    Call state_prop%init(my_lattice,Hams,io_MC%cone)
    dumy=get_charge(my_lattice,Q_neigh)
    qeulerp=dumy(1)
    qeulerm=-dumy(2)
    write(io_status,'(/a)') 'Initialization done'
    write(io_status,'((a,3x,f14.5,a,3x,f9.5))') ' Total energy ', state_prop%E_total, 'eV    Topological charge ', (qeulerp+qeulerm)/pi*0.25d0
    write(io_status,'(a/)') '---------------------'
    
    iT_global=[1,io_MC%n_Tsteps]
    Call distrib_int_index(com,iT_global,iT_local,displ,cnt)
    Call measure_init(measure,it_local,iT_global,com) 
    i_measure=0
    Do i_kT=iT_local(1),iT_local(2)
        i_measure=i_measure+1
        kt=kt_all(i_kT)
        measure(i_measure)%kt=kt
        call Relaxation(my_lattice,io_MC,N_cell,state_prop,qeulerp,qeulerm,kt,hams,Q_neigh)
        !Monte Carlo steps, calculate the values
        do i_MC=1+io_MC%restart_MC_steps,io_MC%Total_MC_Steps+io_MC%restart_MC_steps
            !Monte Carlo steps for independency
            Do i_relax=1,io_MC%T_auto*N_cell
                Call MCstep(my_lattice,io_MC,N_cell,state_prop,kt,Hams)
            End do
            Call MCstep(my_lattice,io_MC,N_cell,state_prop,kt,Hams) ! at least one step
            Call measure_add(measure(i_measure),my_lattice,state_prop,Q_neigh,flat_nei,indexNN) 
        end do ! over i_MC
   
        Call measure_eval(measure(i_measure),io_MC%Cor_log,N_cell)
        if(io_simu%io_Xstruct) Call write_config('MC',kt,my_lattice,filen_kt_acc)
        Write(io_status,'(I6,a,I6,a,f8.4,a,/)')  i_kT, ' nd step out of ', io_MC%n_Tsteps,' steps. T=', kT/k_B,' Kelvin'
    end do !over i_kT

    Call measure_gather(measure,com,displ,cnt)
    if (io_simu%io_spstmL) call spstm  ! The program of Tobias is used only at last iteration
    Call measure_print_thermo(measure,com) ! write EM.dat output
end subroutine montecarlo
end module
