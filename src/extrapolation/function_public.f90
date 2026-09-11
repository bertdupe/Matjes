module m_function_public
use m_function_base
use m_func_linear
use m_func_linear_2d
use m_io_files_utils
use m_io_utils
implicit none

public

contains

subroutine choose_func(func_out)
  class(function_base),allocatable,intent(out) :: func_out

  integer :: io_param
  character(len=50) :: function_name=""

  io_param=open_file_read('input')
  call get_parameter(io_param,'input','function_fit',function_name)
  call close_file('input',io_param)

  if (len_trim(function_name).eq.0) return

  select case (trim(function_name))
     case('linear')
        allocate(func_linear::func_out)
        call func_out%set(2,'linear')
     case('linear_2d')
        allocate(func_linear_2d::func_out)
        call func_out%set(3,'linear_2d')
     case default
        stop 'function not coded'
  end select

end subroutine



end module
