module m_parallel_tempering
use m_paratemp
implicit none
contains
subroutine parallel_tempering(my_lattice,io_simu,ext_param,Hams,com)
    use mpi_distrib_v
    use mpi_basic
    use m_H_public
    use m_topocharge_all
    use m_set_temp
    use m_average_MC
    use m_constants, only : k_b
    use m_vector, only : norm
    use m_store_relaxation
    use m_check_restart
    use m_createspinfile
    use m_derived_types, only : t_cell,io_parameter,simulation_parameters
    use m_derived_types, only : lattice
    use m_topo_commons, only : neighbor_Q,get_charge
    use m_convert
    use m_io_files_utils
    use m_MC_io
    use m_MCstep
    use m_mc_exp_val
    use m_mc_track_val,only: track_val
    use m_fluct, only: fluct_parameters,init_fluct_parameter
    use m_write_config
    use m_io_files_utils, only: close_file

    type(lattice), intent(inout) :: my_lattice
    type(io_parameter), intent(in) :: io_simu
    type(simulation_parameters), intent(in) :: ext_param
    class(t_H), intent(in)                  :: Hams(:)
    class(mpi_type),intent(in)              :: com
    ! internal variable
    type(MC_input)             :: io_MC
    ! slope of temperature sets
    integer :: j_optset
    ! slope of the MC
    integer :: i_MC,i_relax
    !size parameters 
    integer :: N_adjT  !number of temperature adjust cycles (outermost loop)
    integer :: NT_global ! number of global temperatures 
    integer :: NT_local  ! number of local temperatures 
    integer :: size_collect !required local size for parameters that are collected through mpi

    ! internal variable
    ! keep the information of the energies of each replicas before each swap
    real(8), allocatable :: E_temp(:)
    ! considered temperature
    real(8) :: kT,kTfin,kTini
    ! autocorrelation for the parallel tempering and number of relaxation steps
    integer :: autocor_steps,relaxation_steps,total_relax_steps
    ! slope for the number of images and temperatures
    integer :: i_temp
    ! dummy slopes
    integer :: i,io_EM,N_cell,n_sizerelax
    logical :: i_optTset

    !mpi parameters
    type(mpi_type)              :: com_inner    !communicator for inner parallelization, so far ignored (get_two_level_comm) has to be checked first)
    type(mpi_distv)             :: com_outer    !communicator for parallelization of Temperatures
    !integer                     :: iT_global(2),iT_local(2)
    !integer                     :: cnt(com%Np)
    !integer                     :: displ(com%Np)

    type(track_val),allocatable :: state_prop(:)
    type(exp_val),allocatable   :: measure(:)
    
    type(lattice),allocatable,dimension(:)     :: lattices
    type(fluct_parameters)  :: fluct_val
    integer,allocatable     :: Q_neigh(:,:)
    type(paratemp_track)    :: para_track


    
    !input parameters that most probably need some external input, which was not supplied in the version this modification is based on
    N_adjT=1 
    i_optTset=.False.

    !set some input parameters
    call rw_MC(io_MC)
    NT_global=io_MC%n_Tsteps
    n_sizerelax=io_MC%n_sizerelax
    autocor_steps=io_MC%T_auto
    N_cell=my_lattice%Ncell
    kTini=ext_param%ktini
    kTfin=ext_param%ktfin
    ! relaxation variables
    relaxation_steps=n_sizerelax
    total_relax_steps=relaxation_steps*(2**N_adjT-1)  !N_adjTS UNITIALIZED ALSO IN OLD VERSION???
    if(com%ismas)then
        write(6,'(/,a,I6,a,/)') "setup for parallel tempering"
        write(6,'(a,I6,a,/)') "you are calculating",NT_global," temperatures"
        write(6,'(/,a)') 'parallel tempering selected'
        write(6,'(a,I10,a,I10,/)') 'the number of relaxation steps will go from', relaxation_steps, ' to ', total_relax_steps
    endif

    !THIS SHOULD ONLY BE CALCULATED ON ONE MPI THREAD, OTHERWISE THE OUTPUT IS MADNESS
    !get some additional data arrays necessary on all threads for the calculation
    Call neighbor_Q(my_lattice,Q_neigh)
    Call init_fluct_parameter(fluct_val,my_lattice,io_MC%do_fluct)

    !find out how to distribute the threads
!    IT_global=[1,io_MC%n_Tsteps]
!    Call distrib_int_index(com,iT_global,iT_local,displ,cnt)
    Call get_two_level_comm(com,NT_global,com_outer,com_inner) !will not work for inner communicator so far

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
        call ini_temp(para_track%kt_all,kTfin,kTini,NT_global,io_simu%io_warning)
        do i=1,NT_global
            measure(i)%kt=para_track%kt_all(i)
        enddo
        ! open the data file where the thermodynamical quantities should be written
        Call measure_print_thermo_init(io_EM)
    endif

    !Scatter the measure-tasks to all treads
    Call scatterv(measure,com_outer)

    !initialize all local values
    do i_temp=1,NT_local
        Call state_prop(i_temp)%init(lattices(i_temp),Hams,io_MC%cone)
    enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   START THE ACTUAL CALCULATION LOOP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    do j_optset=1,N_adjT
        !  T_relax is T_relax_temp, this is because
        !  one might want to thermalize the first iteration more as compared with the j_optset>2
        if (j_optset.gt.1) autocor_steps=io_MC%T_relax
    
        ! Relaxation of the System
        do i_relax=1,relaxation_steps
            !do loop over local temperatures
            do i_temp=1,NT_local
                kt=measure(i_temp)%kt

                !does it make more sense to do big relaxation step at the very beginning?

                Do i_MC=1,autocor_steps*N_cell
                    Call MCstep(lattices(i_temp),io_MC,N_cell,state_prop(i_temp),kt,Hams)
                enddo
                Call MCstep(lattices(i_temp),io_MC,N_cell,state_prop(i_temp),kt,Hams)! In case autocor_steps set to zero at least one MCstep is done

                Call measure_add(measure(i_temp),lattices(i_temp),state_prop(i_temp),Q_neigh,fluct_val) !add up relevant observable
                E_temp(i_temp)=state_prop(i_temp)%E_total !save the total energy of each replicas
            enddo

            ! reorder the temperature and the replicas
            ! The replicas are actually attached to the processors so they stay remote and can not be rearranged
            ! this part makes sure that the temperature are on the processors of the replicas and that Image_temp points to the right temperature
            ! and then scatter the data on the good remote processors
            Call reorder_temperature(E_temp,measure,com_outer,i_relax,para_track)
        enddo
    
        if (io_simu%io_Xstruct) then !write out config
            do i_temp=1,NT_local
                Call write_config('paratemp_',measure(i_temp)%kt/k_b,lattices(i_temp))!,filen_kt_acc)
            enddo
        endif
    
        !at the end of the relaxation and the image loops, calculate the thermodynamical quantitites
        if(com%ismas)then
            do i=1,NT_global
                Call measure_eval(measure(i),io_MC%Cor_log,N_cell)
            enddo
            Call measure_print_thermo(measure,com,io_EM)
            write(io_EM,'(/)')
            !THIS ORIGINALLY CALLED IRRESPECTIVE OF i_optTset, PUT IT BACK IN IF YOU KNOW WHAT THIS DOES
            !I ONLY MADE A WRAPPER USING PARA_TRACK, EVERYTHING ELSE IS UNCHANGED IN THERE
            if(i_optTset) call calculate_diffusion(para_track,relaxation_steps,i_optTset)! calculate the current, diffusivity... for the paratemp
        endif

        if ((i_optTset).and.(j_optset.ne.N_adjT)) then ! update the temperature set
            if(com%ismas)then
                para_track%kT_all=para_track%kt_updated
                do i_temp=1,NT_global
                    measure(para_track%image_temp(i_temp))%kt=para_track%kt_all(i_temp)
                enddo
                write(6,'(/,a,I5,a,I5)') '#-------- parallel tempering step', j_optset ,' of ', N_adjT
                write(6,'(a,I10,a,I10,/)') '#-------- number of relaxation steps goes from', relaxation_steps ,' --- to --- ', relaxation_steps*2
            endif
            relaxation_steps=relaxation_steps*2
            Call scatterv(measure,com_outer)    !update remote temperatures
        endif
    enddo !over j_MC->N_adjT
    if (com%ismas)then
        write(6,'(a)') 'parallel tempering is done'
        call close_file('EM.dat',io_EM)
    endif
end subroutine parallel_tempering
end module
