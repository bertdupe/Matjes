module m_hessian

contains

subroutine hessian(lat,io_simu,ext_param,H,H_res)
    use m_topo_commons, only: get_charge, neighbor_Q
    use m_update_time, only: update_time,init_update_time
    use m_constants, only : pi,k_b
    use m_energyfield, only : write_energy_field
    use m_createspinfile, only: CreateSpinFile
    use m_user_info, only: user_info
    use m_excitations
    use m_solver_commun, only: get_integrator_field, get_propagator_field,select_propagator
    use m_topo_sd, only: get_charge_map
    use m_solver_order,only: get_dt_mode
    use m_tracker, only: init_tracking,plot_tracking
    use m_print_Beff, only: print_Beff
    use m_precision, only: truncate
    use m_write_config, only: write_config
    use m_energy_output_contribution, only:Eout_contrib_init, Eout_contrib_write
    use m_solver_order,only : get_Dmode_int
    use,intrinsic :: iso_fortran_env, only : output_unit, error_unit
    use m_torque_measurements
!$  use omp_lib

    ! input
    type(lattice), intent(inout)    :: lat
!    type(dyna_input),intent(in)     :: io_dyn
    type(io_parameter), intent(in)  :: io_simu
    type(simulation_parameters), intent(in) :: ext_param
    type(hamiltonian),intent(inout)     ::  H,H_res
!    type(mpi_type),intent(in)       :: comm

! internal variables
    real(8) :: Edy
    integer :: N_cell,dim_m,dim_u
    real(8),allocatable,dimension(:),target :: Beff(:)
    real(8),pointer,contiguous              :: Beff_v(:,:)
    real(8),pointer,contiguous              :: Beff_3(:,:)

    dim_m=lat%M%dim_mode
    dim_u=lat%U%dim_mode
    N_cell=lat%Ncell

    allocate(Beff(dim_mode*N_cell),source=0.0d0)
    Beff_v(1:dim_mode,1:N_cell)=>Beff
    Beff_v(1:dim_mode,1:N_cell)=>Beff
    Beff_3(1:3,1:N_cell*(dim_mode/3))=>Beff

    Edy=H%energy(mag_lattice)
    write(output_unit,'(a,2x,E20.12E3)') 'Initial Total Energy (eV/f.u.)',Edy/real(N_cell,8)

    Call H%get_eff_field(lat_1,Beff,1)

end subroutine

end module m_hessian
