module m_exc_t_heaviside
!this contains the heaviside function for the excitation shape
!it is not really a heaviside function, but the name was chosen like that and I won't change it at this point
use m_exc_t_base, only: excitation_t
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit
implicit none
private
public :: heaviside_read_string

contains

subroutine heaviside_read_string(this,str)
    class(excitation_t),intent(inout)   :: this
    character(len=*),intent(in)             :: str

    integer             :: stat
    character(len=100)  :: dummy_name
    character(len=30)   :: shape_name

    read(str,*,iostat=stat)  dummy_name, shape_name, this%t_start, this%t_end, this%dim_mode
    if(stat/=0)then
        write(error_unit,'(A)') "ERROR, Failed to read excitation shape from line:"
        write(error_unit,'(A)') str
        write(error_unit,'(A)') "The error happens when reading the start end and time (2 reals) and the dim_mode (integer), after the first 2 strings"
        STOP
    endif
    allocate(this%real_var(this%dim_mode*2))
    read(str,*,iostat=stat)  dummy_name, shape_name, this%t_start, this%t_end, this%dim_mode, this%real_var
    if(stat/=0)then
        write(error_unit,'(A)') "ERROR, Failed to read excitation shape from line:"
        write(error_unit,'(A)') str
        write(error_unit,'(A,I3,A)') "The error happens when reading the ", size(this%real_var)," reals after 2 strings, 2 reals, and 1 int"
        STOP
    endif
    this%get_shape=>shape_heaviside
    this%print_t=>print_t_heaviside
end subroutine

function shape_heaviside(this,time) result(val)
    class (excitation_t),intent(in) :: this
    real(8),intent(in)                  :: time
    real(8)                             :: val(this%dim_mode)

    val=heaviside(time,&
                  this%dim_mode,&
                  this%t_start,&
                  this%t_end,&
                  this%real_var(1              :  this%dim_mode),&
                  this%real_var(this%dim_mode+1:2*this%dim_mode))
end function

function heaviside(time,dim_mode,t_start,t_end,val_start,val_end)result(val)
    real(8), intent(in) :: time,t_start,t_end
    integer,intent(in)  :: dim_mode
    real(8), intent(in) :: val_start(dim_mode),val_end(dim_mode)
    real(8)             :: val(dim_mode)

    val=0.d0
    if ((time>=t_start).and.(time<t_end)) val=val_start
    if (time>=t_end) val=val_end(:)
end function

subroutine print_t_heaviside(this,io)
    class(excitation_t),intent(in)  :: this
    integer,intent(in)              :: io
    character(len=10)   ::  dim_mode

    write(dim_mode,'(I10)') this%dim_mode
    write(io,'(3X,A)') "Time-shape: heaviside (obsolete)"
    write(io,'(6X,A)') "Obsolete shape which has a constant I0 within time range and a constant I1 aftter the time-window"
    write(io,'(6X,A)') "Parameters:"
    write(io,'(9X,A,2(F14.4,A))') "time range:    : [",this%t_start,",",this%t_end," ) fs"
    write(io,'(9X,A,'//dim_mode//'(E16.8))') "const. I0      : ",this%real_var(1:this%dim_mode)
    write(io,'(9X,A,'//dim_mode//'(E16.8))') "const. I1      : ",this%real_var(1+this%dim_mode:this%dim_mode*2)
end subroutine

end module 
