module m_random_public
use m_random_number_library
use m_random_mkl
use m_randist
use m_random_base
use m_derived_types
use m_io_utils
use m_io_files_utils
use, intrinsic :: iso_fortran_env, only : error_unit,output_unit
implicit none

integer :: mode = 1
public

contains

subroutine get_ran_type(ran_out)
    class(ranbase),intent(inout),allocatable  :: ran_out

    call read_ran_type('input')

    select case(mode)
        case(1)
           write(output_unit,'(a/)') " Using the internal random number generator"
           allocate(ranint::ran_out)

        case(2)
#if defined CPP_MKL
           write(output_unit,'(a/)') " Using the MKL random number generator"
           allocate(ranmkl::ran_out)
           call ran_out%read_option()
#else
           ERROR STOP "CANNOT USE mkl random number generator: compilation without mkl (CPP_MKL)"
#endif
!        case(3)
!           write(output_unit,'(a/)') " Using the Julia random number generator"
!           allocate(julia_ran::ran_out)
!           call ran_out%read_option()

        case(-1)
            ERROR STOP "cannot found the chosen random number generator"

        case default
            ERROR STOP "unknown random number generator"
    end select

    call ran_out%init_seed()

end subroutine


subroutine read_ran_type(fname)
character(len=*), intent(in) :: fname

integer :: io_input

io_input=open_file_read(fname)
call get_parameter(io_input,fname,'rnd_num_type',mode)
call close_file(fname,io_input)

end subroutine

end module
