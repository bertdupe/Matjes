module m_hessian
use m_derived_types, only : lattice
use m_derived_types, only : io_parameter,simulation_parameters
use m_hamiltonian_collection, only: hamiltonian
use m_hessian_io, only: hessian_input
use mpi_basic, only: mpi_type

contains

subroutine hessian(lat,io_simu,ext_param,H,H_res)
    use m_topo_commons, only: get_charge, neighbor_Q
    use m_constants, only : pi,k_b
    use m_energyfield, only : write_energy_field
    use m_createspinfile, only: CreateSpinFile
    use m_user_info, only: user_info
    use m_topo_sd, only: get_charge_map
    use m_print_Beff, only: print_Beff
    use m_precision, only: truncate
    use m_write_config, only: write_config
    use m_energy_output_contribution, only:Eout_contrib_init, Eout_contrib_write
    use m_solver_order,only : get_Dmode_int
    use,intrinsic :: iso_fortran_env, only : output_unit, error_unit
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
    integer :: N_cell,dim_m,dim_u,dim_mode
    real(8),allocatable,dimension(:),target :: Beff(:)
    real(8),pointer,contiguous              :: Beff_v(:,:)
    real(8),pointer,contiguous              :: Beff_3(:,:)

    dim_m=lat%M%dim_mode
    dim_u=lat%U%dim_mode
    dim_mode=dim_m+dim_u
    N_cell=lat%Ncell

    allocate(Beff(dim_mode*N_cell),source=0.0d0)
    Beff_v(1:dim_mode,1:N_cell)=>Beff
    Beff_v(1:dim_mode,1:N_cell)=>Beff
    Beff_3(1:3,1:N_cell*(dim_mode/3))=>Beff

    Edy=H%energy(lat)
    write(output_unit,'(a,2x,E20.12E3)') 'Initial Total Energy (eV/f.u.)',Edy/real(N_cell,8)

    Call H%get_eff_field(lat,Beff,1)

end subroutine

end module m_hessian
