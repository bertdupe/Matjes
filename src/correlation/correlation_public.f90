module m_correlation_public
use m_correlation_int
use m_correlation_mkl
use m_correlation_base
use m_io_utils
use m_io_files_utils
use, intrinsic :: iso_fortran_env, only : error_unit,output_unit
implicit none

integer :: mode = 1
public

contains

subroutine get_correlation_type(correlation_out)
    class(corre_base),intent(inout),allocatable  :: correlation_out

    call read_correlation_type('input')

    select case(mode)
        case(1)
           write(output_unit,'(a/)') " Using the internal correlation algorithm"
           allocate(corre_int::correlation_out)

        case(2)
#if defined CPP_MKL
           write(output_unit,'(a/)') " Using the MKL correlation algorithm"
!           allocate(corre_mkl::correlation_out)
#else
           ERROR STOP "CANNOT USE mkl rcorrelation algorithm: compilation without mkl (CPP_MKL)"
#endif

        case(-1)
            ERROR STOP "cannot find the chosen correlation algorithm"

        case default
            ERROR STOP "unknown correlation algorithm"
    end select

end subroutine


subroutine read_correlation_type(fname)
character(len=*), intent(in) :: fname

integer :: io_input

io_input=open_file_read(fname)
call get_parameter(io_input,fname,'correlation_type',mode)
call close_file(fname,io_input)

end subroutine

end module
