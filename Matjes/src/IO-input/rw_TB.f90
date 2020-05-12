module m_rw_TB


private
public :: rw_TB

contains

subroutine rw_TB()
use m_io_utils
use m_io_files_utils
implicit none
integer :: io_input
real(kind=8) :: t_1,mu

io_input=open_file_read('input')

call get_parameter(io_input,'input','t_1',t_1)
call get_parameter(io_input,'input','mu',mu)
call close_file('input',io_input)

write(*,*) t_1,mu
stop
end subroutine

end module  m_rw_TB
