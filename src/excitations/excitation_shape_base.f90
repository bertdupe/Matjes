module m_excitation_shape_base
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit

private
public  ::  excitation_shape

type excitation_shape
    character(len=30)                   :: name='uninit' !name for identification     
    real(8)                             :: t_start=0.0d0,t_end=100.0d0
    integer                             :: dim_mode=0   !size of the inner dim_mode for the considered parameter (typically 1 or 3 for scalar/vectors)
    integer                             :: Nvar=0
    real(8),allocatable                 :: real_var(:)
    procedure(int_shape),pointer,pass   :: get_shape=>uninit
contains
    procedure   :: unset
!    procedure   ::   read_string
end type

abstract interface
    function int_shape(this,time)result(val)
        import excitation_shape
        class(excitation_shape),intent(in)  :: this
        real(8), intent(in)                 :: time
        real(8)                             :: val(this%dim_mode)
    end function
end interface

contains

subroutine unset(this)
    class(excitation_shape),intent(inout)   :: this

    this%name='uninit'
    this%t_start=0.0d0
    this%t_end=100.0d0
    this%dim_mode=0
    this%Nvar=0
    deallocate(this%real_var)
    this%get_shape=>uninit
end subroutine


function uninit(this,time)result(val)
    class(excitation_shape),intent(in)  :: this
    real(8), intent(in)                 :: time
    real(8)                             :: val(this%dim_mode)

    Error STOP "Tryint to use an uninitialized excitation_shape"
end function
end module
