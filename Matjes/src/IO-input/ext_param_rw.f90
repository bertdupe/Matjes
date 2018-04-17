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
implicit none
type(simulation_parameters), intent(inout) :: ext_param
! internal variables
integer :: io_input

io_input=open_file_read('input')

call close_file('input',io_input)

end subroutine ext_param_rw
