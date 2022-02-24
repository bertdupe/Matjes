module m_montecarlo
use m_H_public
use m_derived_types, only : lattice,io_parameter,simulation_parameters
use mpi_basic
use m_MC_io
use mpi_util
use m_hamiltonian_collection, only: hamiltonian
implicit none
contains
!
! Routine that does the Monte Carlo (and not the parallel tempering)
!

subroutine montecarlo(lat,io_simu,ext_param,H,com_all)
    !wrapper before Montecarlo_run to bcast MPI-stuff and do some further preparations
    type(lattice), intent(inout)                :: lat
    type(io_parameter), intent(inout)           :: io_simu
    type(simulation_parameters), intent(inout)  :: ext_param
    type(hamiltonian),intent(inout)             :: H
    type(mpi_type),intent(in)                   :: com_all

    type(MC_input)      :: io_MC        !input parameters for MonteCarlo calculation

    Call lat%bcast(com_all)
    Call io_simu%bcast(com_all)
    Call ext_param%bcast(com_all)
    Call H%bcast(com_all)

    if(com_all%ismas) call io_MC%read_file()
    Call io_MC%bcast(com_all)

    Call montecarlo_run(lat,io_MC,io_simu,ext_param,H,com_all)
end subroutine

subroutine montecarlo_run(lat,io_MC,io_simu,ext_param,H,com_all)
    use, intrinsic :: iso_fortran_env, only : output_unit
    use mpi_distrib_v
    use m_constants, only : k_b,pi
    use m_topo_commons, only : neighbor_Q,get_charge
    use m_MCstep
    use m_relaxation
    use m_write_config
    use m_set_temp
    use m_mc_exp_val
    use m_mc_track_val,only: track_val
    use m_fluct, only: fluct_parameters,init_fluct_parameter, print_fluct_spatial, eval_fluct_spatial
    use m_mc_therm_val
    use m_work_ham_single, only: work_ham_single

    type(lattice), intent(inout)            :: lat
    type(MC_input), intent(in)              :: io_MC        !input parameters for MonteCarlo calculation
    type(io_parameter), intent(in)          :: io_simu
    type(simulation_parameters), intent(in) :: ext_param
    type(hamiltonian),intent(inout)         :: H
    type(mpi_type),intent(in)               :: com_all
    
    !!!!!!!!!!!!!!!!!!!!!!!
    ! internal variables
    !!!!!!!!!!!!!!!!!!!!!!!
    ! slope of the MC
    integer                     :: i_relax,i_kT,i_MC
    integer                     :: N_spin       !number of spins in the unit-cell
    integer                     :: NT_global    !global number of temperatures calculated
    integer                     :: NT_local     !local(on mpi-thread) number of temperatures calculated
    integer                     :: size_collect !required local size for parameters that are collected through mpi
    integer                     :: Nstep_auto   !number of MCsteps between each evaluation of the state

    type(track_val)             :: state_prop   !type which saves the current state of the system
    type(exp_val),allocatable   :: measure(:)   !type to keep track of temperature & all expectation values for each temperature
    type(therm_val),allocatable :: thermo(:)    !type to store  thermodynamic properties for each temperature
    integer,allocatable         :: Q_neigh(:,:) !neighbors for topological charge calculation
    type(fluct_parameters)      :: fluct_val    !parameters for fluctuation calculation    
    complex(8),allocatable      :: fluct_spatial(:,:)   !save the spatial distribution of the fluctuations
    !mpi parameters
    type(mpi_type)              :: com_inner    !communicator for inner parallelization, so far ignored (get_two_level_comm)
    type(mpi_distv)             :: com_outer    !communicator for parallelization of Temperatures
    !unimportant local parameters
    real(8),allocatable         :: kt_all(:)
    integer,parameter           :: io_status=output_unit
    integer                     :: filen_kt_acc(2)  !some control about file-output name

    type(work_ham_single)       :: work !type containing work arrays for single energy evaluation
    ! initialize the variables
    N_spin=lat%Ncell*lat%nmag
    NT_global=io_MC%n_Tsteps
    Nstep_auto=max(1,io_MC%T_auto*N_spin)   !max for  at least one step 

    if(com_all%ismas) write(io_status,'(/,a,I6,a,I6,a/)') "you are calculating",io_MC%n_Tsteps," temperatures on ", com_all%NP," procs"

    !find out how to parallelize
    Call get_two_level_comm(com_all,NT_global,com_outer,com_inner)
    NT_local=com_outer%cnt(com_outer%id+1)

    !initialize all values that exist on all threads, but need to be gathered at on master occasionally
    size_collect=NT_local
    if(com_all%isMas) size_collect=NT_global
    allocate(measure(size_collect))
    if(com_inner%ismas)then
        allocate(thermo(size_collect))
    endif

    !initialize necessary things on the master thread (temperatures,Q_neighbors)
    if(com_all%ismas)then
        if(io_MC%expval_read)then
            Call measure_read(measure)
        else
            call ini_temp(kt_all,ext_param%ktfin,ext_param%ktini,NT_global,io_simu%io_warning)
            measure%kt=kt_all
            deallocate(kt_all)
        endif
        filen_kt_acc=[max(int(log10(maxval(measure%kt))),1),5]
        if(io_simu%calc_topo) Call neighbor_Q(lat,Q_neigh)!get neighbors to topological charge calculation
    endif

    !bcast parameters to all threads
    Call bcast(filen_kt_acc,com_all)
    Call bcast_alloc(Q_neigh,com_all)
    Call init_fluct_parameter(fluct_val,lat,io_MC%do_fluct,io_MC%fluct_dir)  !check if this can be bcasted as well

    !Scatter the measure-tasks to all inner master threads (this distributes the temperatures) 
    if(com_inner%ismas)  Call measure_scatterv(measure,com_outer) 

    !intialize tracking variables (total energy, magnetization sum)
    Call state_prop%init(lat,H,io_MC) 

    Call H%get_single_work(1,work)  !allocate work arrays for single energy evaluation (1 for magnetism)
    if(io_MC%do_fluct_spatial) allocate(fluct_spatial(fluct_val%get_nneigh(),lat%Ncell),source=cmplx(0.0d0,0.0d0,8))

    Do i_kT=1,NT_local
        Call measure_bcast(measure(i_kt),com_inner)
        lat%T%all_modes=measure(i_kt)%kt !set local temperature field

        if(io_MC%do_fluct_spatial) fluct_spatial=(0.0d0,0.0d0)
        
        !initial relaxation for each temperature to decorrelate from initial state
        call Relaxation(lat,io_MC,N_spin,state_prop,measure(i_kt)%kt,H,Q_neigh,work,logical(com_inner%ismas))

        !Monte Carlo steps, calculate the values
        do i_MC=1,io_MC%Total_MC_Steps
            !Monte Carlo steps for independency
            Do i_relax=1,Nstep_auto
                Call MCstep(lat,io_MC,N_spin,state_prop,measure(i_kt)%kt,H,work)
            End do
            Call measure_add(measure(i_kt),lat,state_prop,Q_neigh,fluct_val) 
            if(io_MC%do_fluct_spatial) Call eval_fluct_spatial(fluct_spatial,lat,fluct_val)
        end do
        Call measure_reduce(measure(i_kt),com_inner)    !combine same T expectation values over inner MPI
        if(com_inner%ismas)then
            Call therm_set(thermo(i_kt),measure(i_kt),io_MC%Cor_log,N_spin)
            if(io_simu%io_Xstruct) Call write_config('MC',measure(i_kt)%kt/k_b,lat,filen_kt_acc)
            Write(io_status,'(I6,a,I6,a,f8.4,a,/)')  i_kT, ' nd step out of ', NT_local,' local steps. T=', measure(i_kt)%kt/k_B,' Kelvin'
#ifdef CPP_DEBUG
            write(*,*) 'energy internal:', state_prop%E_total,' energy total:', H%energy(lat)
#endif      
        endif

        if(io_MC%do_fluct_spatial) call print_fluct_spatial(measure(i_kt)%N_add,measure(i_kt)%kT/k_B,fluct_spatial,com_inner) !write spatial distribution of fluctuations
    end do !over i_kT

    if(com_inner%ismas)then
        if(io_MC%expval_save)then
            Call measure_gatherv(measure,com_outer)
            if(com_all%ismas) Call measure_save(measure)
        endif
        Call thermo_gatherv(thermo,com_outer)
!        if(io_simu%io_spstmL) call spstm  ! The program of Tobias is used only at last iteration (probably not working at all after all these changes)
        if(com_outer%ismas) Call thermo_print(thermo) ! write EM.dat output
    endif
end subroutine montecarlo_run
end module
