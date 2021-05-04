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

    if(com_all%ismas) call rw_MC(io_MC)
    Call bcast(io_MC,com_all)

    Call montecarlo_run(lat,io_MC,io_simu,ext_param,H,com_all)
end subroutine

subroutine montecarlo_run(lat,io_MC,io_simu,ext_param,H,com_all)
    use, intrinsic :: iso_fortran_env, only : output_unit
    use mpi_distrib_v
    use m_constants, only : k_b,pi
!    use m_rw_MC
    use m_topo_commons, only : neighbor_Q,get_charge
    use m_MCstep
    use m_relaxation
    use m_write_config
    use m_set_temp
    use m_mc_exp_val
    use m_mc_track_val,only: track_val
    use m_average_MC, only: get_neighbours
    use m_fluct, only: fluct_parameters,init_fluct_parameter
    use m_mc_therm_val
    use m_get_table_nn,only :get_table_nn

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

    type(track_val)             :: state_prop   !type which saves the current state of the system (
    type(exp_val),allocatable   :: measure(:)   !type to keep track of temperature & all expectation values for each temperature
    type(therm_val),allocatable :: thermo(:)    !type to store  thermodynamic properties for each temperature
    integer,allocatable         :: Q_neigh(:,:)
    type(fluct_parameters)      :: fluct_val    !parameters for fluctuation calculation    
    !mpi parameters
    type(mpi_type)              :: com_inner    !communicator for inner parallelization, so far ignored (get_two_level_comm) has to be checked first)
    type(mpi_distv)             :: com_outer    !communicator for parallelization of Temperatures
    !unimportant local parameters
    real(8),allocatable         :: kt_all(:)
    integer,parameter           :: io_status=output_unit
    integer                     :: filen_kt_acc(2)  !some control about file-output name

    ! initialize the variables
    N_spin=lat%Ncell*lat%nmag
    NT_global=io_MC%n_Tsteps
    if(com_all%ismas) write(io_status,'(/,a,I6,a,/)') "you are calculating",NT_global," temperatures"

    !find out how to parallelize
    Call get_two_level_comm(com_all,NT_global,com_outer,com_inner)
    NT_local=com_outer%cnt(com_outer%id+1)

    !initialize all values that exist on all threads, but need to be gathered at on master occationally
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
        Call neighbor_Q(lat,Q_neigh)!get neighbors to topological charge calculation
    endif

    !bcast parameters to all threads
    Call bcast(filen_kt_acc,com_all)
    Call bcast_alloc(Q_neigh,com_all)
    Call init_fluct_parameter(fluct_val,lat,io_MC%do_fluct)  !check if this can be bcasted as well

    !Scatter the measure-tasks to all inner master threads (this distributes the temperatures) 
    if(com_inner%ismas)  Call measure_scatterv(measure,com_outer) 

    !intialize tracking variables (total energy, magnetization sum)
    !Call state_prop%init(lat,H,io_MC) 
    Call state_prop%init(lat,H,io_MC) 
    
    Do i_kT=1,NT_local
        Call measure_bcast(measure(i_kt),com_inner)
        lat%T%all_modes=measure(i_kt)%kt !set local temperature field

        if(fluct_val%l_use)then
            allocate(measure(i_kt)%MjpMim_ij(fluct_val%get_nneigh(),lat%Ncell),source=cmplx(0.0d0,0.0d0,8))
        endif

        call Relaxation(lat,io_MC,N_spin,state_prop,measure(i_kt)%kt,H,Q_neigh)
        !Monte Carlo steps, calculate the values
        do i_MC=1,io_MC%Total_MC_Steps
            !Monte Carlo steps for independency
            Do i_relax=1,io_MC%T_auto*N_spin
                Call MCstep(lat,io_MC,N_spin,state_prop,measure(i_kt)%kt,H)
            End do
            Call MCstep(lat,io_MC,N_spin,state_prop,measure(i_kt)%kt,H) ! at least one step
            Call measure_add(measure(i_kt),lat,state_prop,Q_neigh,fluct_val) 
        end do ! over i_MC

        Call measure_reduce(measure(i_kt),com_inner)
        if(com_inner%ismas)then
            Call therm_set(thermo(i_kt),measure(i_kt),io_MC%Cor_log,N_spin)
            if(io_simu%io_Xstruct) Call write_config('MC',measure(i_kt)%kt/k_b,lat,filen_kt_acc)
            Write(io_status,'(I6,a,I6,a,f8.4,a,/)')  i_kT, ' nd step out of ', NT_local,' local steps. T=', measure(i_kt)%kt/k_B,' Kelvin'
#ifdef CPP_DEBUG
            write(*,*) 'energy internal:', state_prop%E_total,' energy total:', H%energy(lat)
#endif      
        endif
        if(fluct_val%l_use) call print_spatial_fluct(measure(i_kt),com_all) !write spatial distribution of fluctuations
        if(allocated(measure(i_kt)%MjpMim_ij)) deallocate(measure(i_kt)%MjpMim_ij)
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
