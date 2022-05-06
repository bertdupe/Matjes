module m_exc_t_base
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit

private
public  ::  excitation_t

type excitation_t
    character(len=30)                       :: name='uninit' !name for identification     
    real(8)                                 :: t_start=0.0d0,t_end=100.0d0
    integer                                 :: dim_mode=0   !size of the inner dim_mode for the considered parameter (typically 1 or 3 for scalar/vectors)
    integer                                 :: Nvar=0
    real(8),allocatable                     :: real_var(:)
    procedure(int_shape),pointer,pass       :: get_shape=>uninit
    procedure(int_print_t), pointer,pass    :: print_t=>print_t_uninit
contains
    procedure   :: unset
!    procedure   ::   read_string
end type

abstract interface
    function int_shape(this,time)result(val)
        import excitation_t
        class(excitation_t),intent(in)  :: this
        real(8), intent(in)             :: time
        real(8)                         :: val(this%dim_mode+1)  ! val(1:dim_mode) magnitude
                                                                 ! val(dim_mode+1) complex phase
    end function
    subroutine int_print_t(this,io)
        import excitation_t
        class(excitation_t),intent(in)  :: this
        integer,intent(in)              :: io
    end subroutine 
end interface

contains

subroutine unset(this)
    class(excitation_t),intent(inout)   :: this

    this%name='uninit'
    this%t_start=0.0d0
    this%t_end=100.0d0
    this%dim_mode=0
    this%Nvar=0
    deallocate(this%real_var)
    this%get_shape=>uninit
    this%print_t=>print_t_uninit
end subroutine


function uninit(this,time)result(val)
    class(excitation_t),intent(in)  :: this
    real(8), intent(in)                 :: time
    real(8)                             :: val(this%dim_mode+1)

    val=0.d0
    Error STOP "Tryint to use an uninitialized excitation_t"
end function

subroutine print_t_uninit(this,io)
    class(excitation_t),intent(in)    :: this
    integer,intent(in)                      :: io

    Error STOP "Tryint to print an uninitialized excitation_t"
end subroutine

end module
