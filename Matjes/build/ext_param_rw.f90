!
! Subroutine that reads the external parameters of the Hamiltonian
!
!
!
!

subroutine ext_param_rw(ext_param)
use m_derived_types, only : simulation_parameters
use m_io_utils
use m_io_files_utils
use m_constants, only : k_B,mu_B
implicit none
type(simulation_parameters), intent(inout) :: ext_param
! internal variables
integer :: io_input
real(kind=8) :: kt

io_input=open_file_read('input')

call get_parameter(io_input,'input','H_ext',3,ext_param%H_ext%value)

call get_parameter(io_input,'input','E_ext',3,ext_param%E_ext%value)

call get_parameter(io_input,'input','Tini',kt)
ext_param%ktini%value=kt*k_B
call get_parameter(io_input,'input','Tfin',kt)
ext_param%ktfin%value=kt*k_B

call close_file('input',io_input)

end subroutine ext_param_rw
