module m_exc_t_algebraic
!this contains the algebraic-decay function for the excitation shape
! the equation is:
! val= I0 /(t-t_0)**k
use m_exc_t_base, only: excitation_t
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit
implicit none
private
public :: algebraic_read_string

contains

subroutine algebraic_read_string(this,str)
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
        write(error_unit,'(AI3,A)') "The error happens when reading the ", size(this%real_var)," reals after 2 strings, 2 reals and 1 int"
        STOP
    endif

    associate( I0     => this%real_var(1                :  this%dim_mode  ),&
               t0     => this%real_var(1+this%dim_mode  :  this%dim_mode+1),&
               k      => this%real_var(1+this%dim_mode  :  this%dim_mode+2))

        if(t0(1)>this%t_start.and.t0(1)<this%t_end)then
            write(error_unit,'(A)') "WARNING, in algebraic excitation t-shape the t0 is within (t_start,t_end)"
            write(error_unit,'(A)') "This can lead to divergences which is probably not intended"
            write(error_unit,'(3(AE16.8/))') "t_start: ", this%t_start,"t_end:   ", this%t_end, "t0:     ", t0(1)
            write(error_unit,'(A)') "CHECK INPUT"
            ERROR STOP
        endif
    end associate

    this%get_shape=>shape_algebraic
    this%print_t=>print_t_algebraic
end subroutine

function shape_algebraic(this,time) result(val)
    class (excitation_t),intent(in) :: this
    real(8),intent(in)                  :: time
    real(8)                             :: val(this%dim_mode)

    associate( I0     => this%real_var(1                :  this%dim_mode  ),&
               t0     => this%real_var(1+this%dim_mode  :  this%dim_mode+1),&
               k      => this%real_var(2+this%dim_mode  :  this%dim_mode+2))
    val=algebraic(time,&
                  this%dim_mode,&
                  this%t_start,&
                  this%t_end,&
                  I0,t0(1),k(1))
    end associate
end function

function algebraic(time,dim_mode,t_start,t_end,I0,t0,k)result(val)
    real(8), intent(in) :: time,t_start,t_end
    integer,intent(in)  :: dim_mode
    real(8), intent(in) :: I0(dim_mode)
    real(8), intent(in) :: t0,k
    real(8)             :: val(dim_mode)

    val=0.0d0
    if ((time>=t_start).and.(time<t_end))then
        val=I0/(time-t0)**k
    endif
end function

subroutine print_t_algebraic(this,io)
    use m_constants, only: pi
    class(excitation_t),intent(in)  :: this
    integer,intent(in)              :: io
    character(len=10)   ::  dim_mode

    write(dim_mode,'(I10)') this%dim_mode
    associate( I0     => this%real_var(1                :  this%dim_mode  ),&
               t0     => this%real_var(1+this%dim_mode  :  this%dim_mode+1),&
               k      => this%real_var(2+this%dim_mode  :  this%dim_mode+2))
    write(io,'(3X,A)') "Time-shape: algebraic"
    write(io,'(6X,A)') "Algebraic increase/decay"
    write(io,'(6X,A)') "Equation: I0/(t-t0)**k"
    write(io,'(6X,A)') "Parameters:"
    write(io,'(9X,A,2(F14.4,A))')            "time range:    : [",this%t_start,",",this%t_end," ) fs"
    write(io,'(9X,A,'//dim_mode//'(E16.8))') "amplitude (I0) : ", I0
    write(io,'(9X,A,(F16.4),A)')             "center    (t0) : ", t0(1) , " in fs"
    write(io,'(9X,A,(F16.4))')               "decay      (k) : ", k(1) 
    end associate
end subroutine

end module 
