module m_random_base
use m_io_files_utils
use m_convert
use, intrinsic :: iso_fortran_env, only : error_unit,output_unit
implicit none

private
public :: ranbase

type, abstract :: ranbase
    real(8),allocatable   :: x(:)      ! variable containing all the random numbers
    integer               :: N         ! number of random numbers to be generated
    logical               :: print = .false.   ! print out the numbers
contains
    ! defered type
    procedure(int_seed),deferred          :: init_seed
    procedure(int_getx),deferred          :: get_list
    procedure(int_destroy),deferred       :: destroy
    procedure(int_gextra_list),deferred   :: get_extract_list
    procedure(int_rw_option),deferred     :: read_option
    ! base function
    procedure, NON_OVERRIDABLE     :: init_base
    procedure, NON_OVERRIDABLE     :: print_base
    procedure, NON_OVERRIDABLE     :: extract_list
    procedure, NON_OVERRIDABLE     :: extract_size
end type


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
abstract interface

    subroutine int_seed(this)
       import ranbase
       class(ranbase), intent(in)     :: this
    end subroutine

    subroutine int_getx(this,a,b)
       import ranbase
       class(ranbase), intent(inout)  :: this
       real(8),intent(in)             :: a,b
    end subroutine

    subroutine int_destroy(this)
       import ranbase
       class(ranbase), intent(in)     :: this
    end subroutine

    subroutine int_gextra_list(this,a,b,resu)
       import ranbase
       class(ranbase), intent(inout)  :: this
       real(8),intent(in)             :: a,b
       real(8),intent(inout)          :: resu(:)
    end subroutine

    subroutine int_rw_option(this)
       import ranbase
       class(ranbase), intent(inout)  :: this
    end subroutine

end interface


contains

subroutine extract_list(this,x)
    class(ranbase),intent(in)    :: this
    real(8), intent(inout)       :: x(:)

    integer :: N

    N=size(x)
    x=this%x(1:N)

end subroutine

subroutine extract_size(this,N)
    class(ranbase),intent(in)    :: this
    integer, intent(out)         :: N

    N=this%N

end subroutine

subroutine init_base(this,N)
    class(ranbase),intent(inout) :: this
    integer, intent(in)          :: N

    if (allocated(this%x)) STOP "cannot allocate random number vector - already allocated"

    this%N=N
    allocate(this%x(N),source=0.0d0)

end subroutine


subroutine print_base(this,tag)
    class(ranbase),intent(in) :: this
    integer, intent(in)          :: tag

    integer :: i,io_out
    character(len=100) :: fname

    if (.not.allocated(this%x)) STOP "cannot print random numbers vector - list not allocated"

    fname=convert('rnd_num_',tag)
    io_out=open_file_write(trim(fname))
    do i=1,this%N
       write(io_out,'(E20.12E3)') this%x(i)
    enddo
    call close_file(fname,io_out)

end subroutine

end module