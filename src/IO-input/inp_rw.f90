!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine reads the input file and setup and every variables given in the inp file
!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine inp_rw(io_simu)
use m_constants
use m_derived_types
use m_io_utils
use m_io_files_utils
implicit none
!ccccccccccccccccccccccccccccccccccccccccccc
!In/out variable
type(io_parameter), intent(inout) :: io_simu
! internal variables
integer  :: io_input
! local variables
!ccccccccccccccccccccccccccccccccccccccccccc

io_input=open_file_read('input')

! mpi variables
!call get_parameter(io_input,'input','ghost',i_ghost)
!call get_parameter(io_input,'input','algo_mpi',i_ghost)
!call get_parameter(io_input,'input','nRepProc',nRepProc)

! io variables
call get_parameter(io_input,'input','gra_fft',io_simu%io_fft_Xstruct)
call get_parameter(io_input,'input','gra_topo',io_simu%io_topo)
call get_parameter(io_input,'input','warnings',io_simu%io_warning)
call get_parameter(io_input,'input','gra_log',io_simu%io_Xstruct)
call get_parameter(io_input,'input','gra_StochF',io_simu%io_Tfield)
call get_parameter(io_input,'input','gra_freq',io_simu%io_frequency)
call get_parameter(io_input,'input','gra_Beff',io_simu%io_Beff)
call get_parameter(io_input,'input','Energy_Distrib',io_simu%io_Energy_Distrib)
call get_parameter(io_input,'input','energy_detail',io_simu%io_Energy_detail)
call get_parameter(io_input,'input','print_Econt',io_simu%io_energy_cont)
call get_parameter(io_input,'input','Angle_Distrib',io_simu%io_Angle_Distrib)
call get_parameter(io_input,'input','Field_Distrib',io_simu%io_Field_Distrib)
call get_parameter(io_input,'input','Forces',io_simu%io_Force)
call get_parameter(io_input,'input','tracker',io_simu%io_tracker)
call get_parameter(io_input,'input','calc_topo',io_simu%calc_topo)

call close_file('input',io_input)

end subroutine inp_rw

