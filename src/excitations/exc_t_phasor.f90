module m_exc_t_phasor
!this contains the phasor function for the excitation shape
! the equation is:
! val= I0 * exp(-i*2pi/T*t+phi*pi)
use m_exc_t_base, only: excitation_t
use m_constants
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit
implicit none
private
public :: phasor_read_string

contains

subroutine phasor_read_string(this,str)
    class(excitation_t),intent(inout)   :: this
    character(len=*),intent(in)             :: str

    integer             :: stat
    character(len=100)  :: dummy_name
    character(len=30)   :: shape_name

    read(str,*,iostat=stat)  dummy_name, shape_name, this%t_start, this%t_end, this%dim_mode
    if(stat/=0)then
        write(error_unit,'(A)') "ERROR, Failed to read excitation shape from line:"
        write(error_unit,'(A)') str
        write(error_unit,'(A)') "The error happens when reading the dim_mode (integer), after the first 2 strings"
        STOP
    endif
    allocate(this%real_var(this%dim_mode+2))
    read(str,*,iostat=stat)  dummy_name, shape_name, this%t_start, this%t_end, this%dim_mode, this%real_var
    if(stat/=0)then
        write(error_unit,'(A)') "ERROR, Failed to read excitation shape from line:"
        write(error_unit,'(A)') str
        write(error_unit,'(A,I3,A)') "The error happens when reading the ", size(this%real_var)," reals after 2 strings, 2 reals and 1 int"
        STOP
    endif

    ! associate( I0     => this%real_var(1                :  this%dim_mode  ),&
    !            period => this%real_var(1+this%dim_mode  :  this%dim_mode+1),&
    !            phi    => this%real_var(2+this%dim_mode  :  this%dim_mode+2))
    ! end associate

    this%get_shape=>shape_phasor
    this%print_t=>print_t_phasor
end subroutine

function shape_phasor(this,time) result(val)
    class (excitation_t),intent(in) :: this
    real(8),intent(in)                  :: time
    real(8)                             :: val(this%dim_mode+1)

    associate( I0     => this%real_var(1                :  this%dim_mode  ),&
               period => this%real_var(1+this%dim_mode),&
               phi    => this%real_var(2+this%dim_mode))
    if(time < this%t_start .or. time > this%t_end) then
        val = 0.0d0
    else
        val(1:this%dim_mode) = I0
        val(this%dim_mode+1) = -2.0d0*pi/period * time + phi * pi
    endif
    end associate
end function

subroutine print_t_phasor(this,io)
    use m_constants, only: pi
    class(excitation_t),intent(in)  :: this
    integer,intent(in)              :: io
    character(len=10)   ::  dim_mode

    write(dim_mode,'(I10)') this%dim_mode
    associate( I0     => this%real_var(1                :  this%dim_mode  ),&
               period => this%real_var(1+this%dim_mode  :  this%dim_mode+1),&
               phi    => this%real_var(2+this%dim_mode  :  this%dim_mode+2))
    write(io,'(3X,A)') "Time-shape: phasor"
    write(io,'(6X,A)') "Complex phasor"
    write(io,'(6X,A)') "Equation: I0 * exp(-i*2pi/T*t+phi*pi)"
    write(io,'(6X,A)') "Parameters:"
    write(io,'(9X,A,2(F14.4,A))')            "time range:    : [",this%t_start,",",this%t_end," ) fs"
    write(io,'(9X,A,'//dim_mode//'(E16.8))') "amplitude (I0) : ", I0
    write(io,'(9X,A,(F16.4),A)')             "period     (T) : ", period(1) , " in fs"
    write(io,'(9X,A,(F16.4))')               "phase    (phi) : ", phi(1)
    end associate
end subroutine

end module
