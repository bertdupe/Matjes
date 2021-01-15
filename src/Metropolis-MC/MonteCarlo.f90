module m_montecarlo
implicit none
contains
!
! Routine that does the Monte Carlo (and not the parallel tempering)
!
!subroutine MonteCarlo(lats,mag_motif,io_simu,gra_topo,ext_param)
subroutine montecarlo(lat,io_simu,ext_param,Hams,com_all)
    use, intrinsic :: iso_fortran_env, only : output_unit
    use mpi_distrib_v
    use mpi_util
    use m_constants, only : k_b,pi
    use m_derived_types, only : lattice,t_cell,io_parameter,simulation_parameters
!    use m_rw_MC
    use m_topo_commons, only : neighbor_Q,get_charge
    use m_MCstep
    use m_H_public
    use m_relaxation
    use m_write_config
    use m_set_temp
    use m_mc_exp_val
    use m_mc_track_val,only: track_val
    use m_average_MC, only: get_neighbours
    use m_fluct, only: fluct_parameters,init_fluct_parameter
    use m_MC_io
    use m_get_table_nn,only :get_table_nn

    type(lattice), intent(inout)            :: lat
    type(io_parameter), intent(in)          :: io_simu
    type(simulation_parameters), intent(in) :: ext_param
    class(t_H), intent(in)                  :: Hams(:)
    type(mpi_type),intent(in)               :: com_all
    
    !!!!!!!!!!!!!!!!!!!!!!!
    ! internal variables
    !!!!!!!!!!!!!!!!!!!!!!!
    ! slope of the MC
    integer                     :: i_relax,i_kT,i_MC
    integer                     :: N_spin       !number of spins in the unit-cell
    integer                     :: NT_global    !global number of temperatures calculated
    integer                     :: NT_local     !local(on mpi-thread) number of temperatures calculated
    ! variable for the temperature
!    real(8)                     :: qeulerp,qeulerm,dumy(5)
    type(track_val)             :: state_prop
    ! contribution of the different energy parts
    type(exp_val),allocatable   :: measure(:)   !type to keep track of temperature & all thermodynamic properties for each temperature
    type(MC_input)              :: io_MC        !input parameters from MonteCarlo calculation
    integer,allocatable         :: Q_neigh(:,:)
    type(fluct_parameters)      :: fluct_val    !parameters for fluctuation calculation    
    !mpi parameters
    type(mpi_type)              :: com_inner    !communicator for inner parallelization, so far ignored (get_two_level_comm) has to be checked first)
    type(mpi_distv)             :: com_outer    !communicator for parallelization of Temperatures
    !unimportant local parameters
    real(8),allocatable         :: kt_all(:)
    integer                     :: ierr
    integer,parameter           :: io_status=output_unit
    integer                     :: filen_kt_acc(2)  !some control about file-output name
    integer,allocatable        ::  indexNN(:)
    integer, allocatable :: tableNN(:,:,:,:,:,:)

    ! initialize the variables
    N_spin=lat%Ncell*lat%nmag
    if(com_all%ismas)then
        call rw_MC(io_MC)
        NT_global=io_MC%n_Tsteps
        write(io_status,'(/,a,I6,a,/)') "you are calculating",NT_global," temperatures"

        !get all temperatures (kt_all)
        call ini_temp(kt_all,ext_param%ktfin,ext_param%ktini,NT_global,io_simu%io_warning)
        filen_kt_acc=5
        filen_kt_acc(1)=max(int(log10(maxval(kt_all))),1)

        !initialize measure-array
        allocate(measure(NT_global))
        measure(:)%kt=kt_all(:)
        deallocate(kt_all)
    	Call get_table_nn(lat,1,indexNN,tableNN) !need indexNN to allocate fluctuation arrays in measure

        !get neighbors to topological charge calculation
        Call neighbor_Q(lat,Q_neigh)
    endif

    !find out how to parallelize
    Call bcast(NT_global,com_all)
    Call get_two_level_comm(com_all,NT_global,com_outer,com_inner)
    NT_local=com_outer%cnt(com_outer%id+1)

    !distribute to all threads
    Call bcast(io_MC,com_all)
    Call bcast(filen_kt_acc,com_all)
    Call bcast_alloc(Q_neigh,com_all)
    Call init_fluct_parameter(fluct_val,lat,io_MC%do_fluct)  !check if this can be bcasted as well

    !Scatter the measure-tasks to all treads    
    if(.not.allocated(measure)) allocate(measure(NT_local))
    Call scatterv(measure,com_outer)
    
    !intialize tracking variables (total energy, magnetization sum)
    Call state_prop%init(lat,Hams,io_MC%cone)
    
    Do i_kT=1,NT_local
        lat%T%all_modes=measure(i_kt)%kt !set local temperature field

		allocate(measure(i_kt)%MjpMim_sum(sum(indexNN)))
		allocate(measure(i_kt)%MjpMim_av(sum(indexNN)))
		allocate(measure(i_kt)%MjpMim_ij_sum(sum(indexNN),lat%Ncell))
		allocate(measure(i_kt)%MjpMim_ij_av(sum(indexNN),lat%Ncell))

		measure(i_kt)%MjpMim_sum(:)=cmplx(0.0d0,0.0d0,8)
		measure(i_kt)%MjpMim_av(:)=cmplx(0.0d0,0.0d0,8)
		measure(i_kt)%MjpMim_ij_sum(:,:)=cmplx(0.0d0,0.0d0,8)
		measure(i_kt)%MjpMim_ij_av(:,:)=cmplx(0.0d0,0.0d0,8)

        call Relaxation(lat,io_MC,N_spin,state_prop,measure(i_kt)%kt,hams,Q_neigh)
        !Monte Carlo steps, calculate the values
        do i_MC=1+io_MC%restart_MC_steps,io_MC%Total_MC_Steps+io_MC%restart_MC_steps
            !Monte Carlo steps for independency
            Do i_relax=1,io_MC%T_auto*N_spin
                Call MCstep(lat,io_MC,N_spin,state_prop,measure(i_kt)%kt,Hams)
            End do
            Call MCstep(lat,io_MC,N_spin,state_prop,measure(i_kt)%kt,Hams) ! at least one step
            Call measure_add(measure(i_kt),lat,state_prop,Q_neigh,fluct_val) 
        end do ! over i_MC
   
        Call measure_eval(measure(i_kt),io_MC%Cor_log,N_spin)
        if(io_simu%io_Xstruct) Call write_config('MC',measure(i_kt)%kt/k_b,lat,filen_kt_acc)
		if(fluct_val%l_use) call  print_spatial_fluct(measure(i_kt),com_all) !write spatial distribution of fluctuations
        Write(io_status,'(I6,a,I6,a,f8.4,a,/)')  i_kT, ' nd step out of ', size(measure),' steps. T=', measure(i_kt)%kt/k_B,' Kelvin'

		deallocate(measure(i_kt)%MjpMim_ij_sum)
		deallocate(measure(i_kt)%MjpMim_ij_av)

    end do !over i_kT

    Call gatherv(measure,com_outer)
    if (io_simu%io_spstmL) call spstm  ! The program of Tobias is used only at last iteration
    Call measure_print_thermo(measure,com_all) ! write EM.dat output

end subroutine montecarlo
end module
