module m_fit_public
use m_fit_base
use m_fit_linear
use m_fit_nonlinear
use m_function_public
use m_io_files_utils
use m_io_utils
use, intrinsic :: iso_fortran_env, only : output_unit
implicit none

private
public :: choose_fit

contains

subroutine choose_fit(my_fit)
  class(fit),allocatable,intent(out) :: my_fit

  integer :: io_param
  character(len=50) :: fit_name=""

  io_param=open_file_read('input')
  call get_parameter(io_param,'input','fit_type',fit_name)
  call close_file('input',io_param)

  if (len_trim(fit_name).eq.0) return

  select case (trim(fit_name))
     case('internal_linear')
        allocate(fit_linear::my_fit)

        write(output_unit,'(a)') 'internal linear fit routines selected'

     case('internal_non_linear')
        allocate(fit_nonlinear::my_fit)

        write(output_unit,'(a)') 'internal non-linear fit routines selected'

!     case('mkl')
!        allocate(fit_mkl::my_fit)
     case default
        stop 'fit not coded'
  end select

  call choose_func(my_fit%func)
  call my_fit%func%show(tag=-1.0d0)

end subroutine

end module
