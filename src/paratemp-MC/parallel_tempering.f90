module m_parallel_tempering
use m_paratemp
use mpi_basic
use mpi_util
use m_hamiltonian_collection, only: hamiltonian
use m_derived_types, only : lattice, io_parameter, simulation_parameters
implicit none
contains
subroutine parallel_tempering(my_lattice,io_simu,ext_param,H,com)
    type(lattice), intent(inout)                :: my_lattice
    type(io_parameter), intent(inout)           :: io_simu
    type(simulation_parameters), intent(inout)  :: ext_param
    type(hamiltonian), intent(inout)            :: H
    class(mpi_type), intent(in)                 :: com

    Call my_lattice%bcast(com)
    Call io_simu%bcast(com)
    Call ext_param%bcast(com)
    Call H%bcast(com)

    Call parallel_tempering_run(my_lattice,io_simu,ext_param,H,com)
end subroutine

subroutine parallel_tempering_run(my_lattice,io_simu,ext_param,H,com)
    use mpi_distrib_v
    use m_H_public
    use m_topocharge_all
    use m_set_temp
    use m_average_MC
    use m_constants, only : k_b
    use m_vector, only : norm
    use m_store_relaxation
    use m_check_restart
    use m_createspinfile
    use m_topo_commons, only : neighbor_Q,get_charge
    use m_convert
    use m_io_files_utils
    use m_MC_io
    use m_MCstep
    use m_mc_track_val,only: track_val
    use m_mc_exp_val
    use m_mc_therm_val
    use m_fluct, only: fluct_parameters,init_fluct_parameter
    use m_write_config
    use m_io_files_utils, only: close_file
    use m_work_ham_single, only: work_ham_single

    type(lattice), intent(inout) :: my_lattice
    type(io_parameter), intent(in) :: io_simu
    type(simulation_parameters), intent(in) :: ext_param
    type(hamiltonian),intent(inout)         :: H
    class(mpi_type),intent(in)              :: com
    ! internal variable
    type(MC_input)             :: io_MC
    ! slope of temperature sets
    integer :: j_optset
    ! slope of the MC
    integer :: i_swapT
    integer :: i_MC
    !size parameters 
    integer :: N_adjT  !number of temperature adjust cycles (outermost loop)
    integer :: NT_global ! number of global temperatures 
    integer :: NT_local  ! number of local temperatures 
    integer :: size_collect !required local size for parameters that are collected through mpi

    ! internal variable
    ! keep the information of the energies of each replicas before each swap
    real(8), allocatable :: E_temp(:)
    ! considered temperature
    real(8) :: kT
    ! autocorrelation for the parallel tempering and number of relaxation steps
    integer :: autocor_steps
    integer :: n_swapT  !number of times T might be swapped within T-set
    ! slope for the number of images and temperatures
    integer :: i_temp
    ! dummy slopes
    integer :: i,io_EM,N_spin,n_sizerelax
    logical :: i_optTset

    !mpi parameters
    type(mpi_type)              :: com_inner    !communicator for inner parallelization, so far ignored (get_two_level_comm) has to be checked first)
    type(mpi_distv)             :: com_outer    !communicator for parallelization of Temperatures

    type(track_val),allocatable :: state_prop(:)
    type(exp_val),allocatable   :: measure(:)   !type to keep track of temperature & all expectation values for each temperature
    type(therm_val),allocatable :: thermo(:)    !type to store  thermodynamic properties for each temperature
    
    type(lattice),allocatable,dimension(:)     :: lattices
    type(fluct_parameters)  :: fluct_val
    integer,allocatable     :: Q_neigh(:,:)
    type(paratemp_track)    :: para_track

    type(work_ham_single)       :: work !type containing work arrays for single energy evaluation
    
    !input parameters that most probably need some external input, which was not supplied in the version this modification is based on
    i_optTset=.True.

    !set some input parameters
    if(com%ismas) call rw_MC(io_MC)
    Call io_MC%bcast(com)
    N_adjT=io_MC%N_Topt
    NT_global=io_MC%n_Tsteps
    n_sizerelax=io_MC%n_sizerelax
    autocor_steps=io_MC%T_auto
    N_spin=my_lattice%Ncell*my_lattice%nmag
    n_swapT=n_sizerelax
    if(com%ismas)then
        write(6,'(/,a,I6,a,/)') "setup for parallel tempering"
        write(6,'(a,I6,a,/)') "you are calculating",NT_global," temperatures"
        write(6,'(/,a)') 'parallel tempering selected'
        write(6,'(a,I10,a,I10,/)') 'the number of relaxation steps will go from', n_swapT, ' to ', n_swapT*(2**N_adjT-1)
    endif


    !find out how to distribute the threads
    Call get_two_level_comm(com,NT_global,com_outer,com_inner) !will not work for inner communicator so far
    if(com_outer%NP>NT_global) STOP "should not use more mpi-threads that parallel tempering temperatures"

    !initialize all values that will only exist locally on each MPI-thread
    NT_local=com_outer%cnt(com%id+1)
    allocate(lattices(NT_local))
    allocate(state_prop(NT_local))
    do i_temp=1,NT_local
        Call my_lattice%copy(lattices(i_temp))
    enddo

    !initialize all values that exist on all threads, but need to be gathered at on master occationally
    size_collect=NT_local
    if(com%isMas) size_collect=NT_global
    allocate(E_temp(size_collect),source=0.0d0)
    allocate(measure(size_collect))

    !initialize necessary things on the master thread (temperatures,io-stuff)
    if(com%ismas)then
        Call alloc_paratemp_track(para_track,NT_global)
        call ini_temp(para_track%kt_all,ext_param%ktfin,ext_param%ktini,NT_global,io_simu%io_warning)
        do i=1,NT_global
            measure(i)%kt=para_track%kt_all(i)
        enddo
        ! open the data file where the thermodynamical quantities should be written
        Call thermo_print_init(io_EM)
        allocate(thermo(NT_global))
        Call neighbor_Q(my_lattice,Q_neigh)
    endif
   
    !broadcast and calculate locally necessary values for further calculations
    Call bcast_alloc(Q_neigh,com)
    Call init_fluct_parameter(fluct_val,my_lattice,io_MC%do_fluct,io_MC%fluct_dir)  !check if this can be bcasted as well

    !Scatter the measure-tasks to all treads
    Call measure_scatterv(measure,com_outer)

    !initialize all local values
    do i_temp=1,NT_local
        Call state_prop(i_temp)%init(lattices(i_temp),H,io_MC)
    enddo

    Call H%get_single_work(1,work)  !allocate work arrays for single energy evaluation (1 for magnetism)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   START THE ACTUAL CALCULATION LOOP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    do j_optset=1,N_adjT
        !  T_relax is T_relax_temp, this is because
        !  one might want to thermalize the first iteration more as compared with the j_optset>2
        if (j_optset.gt.1) autocor_steps=io_MC%T_relax
    
        ! Relaxation of the System
        do i_swapT=1,n_swapT
            !do loop over local temperatures
            do i_temp=1,NT_local
                kt=measure(i_temp)%kt

                !does it make more sense to do big relaxation step at the very beginning?
                Do i_MC=1,autocor_steps*N_spin
                    Call MCstep(lattices(i_temp),io_MC,N_spin,state_prop(i_temp),kt,H,work)
                enddo
                Call MCstep(lattices(i_temp),io_MC,N_spin,state_prop(i_temp),kt,H,work)! In case autocor_steps set to zero at least one MCstep is done

                Call measure_add(measure(i_temp),lattices(i_temp),state_prop(i_temp),Q_neigh,fluct_val) !add up relevant observable
                E_temp(i_temp)=state_prop(i_temp)%E_total !save the total energy of each replicas
            enddo

            ! reorder the temperature and the replicas
            ! The replicas are actually attached to the processors so they stay remote and can not be rearranged
            ! this part makes sure that the temperature are on the processors of the replicas and that Image_temp points to the right temperature
            ! and then scatter the data on the good remote processors
            Call reorder_temperature(E_temp,measure,com_outer,i_swapT,para_track)
        enddo
    
        if (io_simu%io_Xstruct) then !write out config
            do i_temp=1,NT_local
                Call write_config('paratemp_',measure(i_temp)%kt/k_b,lattices(i_temp))!,filen_kt_acc)
            enddo
        endif
    
        !at the end of the relaxation and the image loops, calculate the thermodynamical quantitites
        if(com%ismas)then
            do i=1,NT_global
                Call therm_set(thermo(i),measure(i),io_MC%Cor_log,N_spin)
            enddo
            Call thermo_print(thermo,io_EM) ! write EM.dat output
            write(io_EM,'(/)')
            !THIS ORIGINALLY CALLED IRRESPECTIVE OF i_optTset, PUT IT BACK IN IF YOU KNOW WHAT THIS DOES
            !I ONLY MADE A WRAPPER USING PARA_TRACK, EVERYTHING ELSE IS UNCHANGED IN THERE
            if(i_optTset) call calculate_diffusion(para_track,n_swapT,i_optTset)! calculate the current, diffusivity... for the paratemp
        endif

        if ((i_optTset).and.(j_optset.ne.N_adjT)) then ! update the temperature set
            if(com%ismas)then
                para_track%kT_all=para_track%kt_updated
                measure=exp_val()   !reset measure 
                measure(para_track%image_temp)%kt=para_track%kt_all !set new temperatures
                write(6,'(/,a,I5,a,I5)') '#-------- parallel tempering step', j_optset ,' of ', N_adjT
                write(6,'(a,I10,a,I10,/)') '#-------- number of relaxation steps goes from', n_swapT ,' --- to --- ', n_swapT*2
            endif
            n_swapT=n_swapT*2
            Call measure_scatterv(measure,com_outer)    !update remote temperatures
        endif
    enddo !over j_MC->N_adjT
    if (com%ismas)then
        write(6,'(a)') 'parallel tempering is done'
        call close_file('EM.dat',io_EM)
    endif
end subroutine parallel_tempering_run
end module
