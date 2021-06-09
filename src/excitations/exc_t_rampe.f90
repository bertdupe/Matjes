module m_exc_t_rampe
use m_convert
use m_exc_t_base, only: excitation_t
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit

private
public :: rampe_read_string

contains


subroutine rampe_read_string(this,str)
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
        write(error_unit,'(A)') "The error happens when reading the ", size(this%real_var)," reals after 2 strings, 2 reals, and 1 int"
        STOP
    endif
    this%get_shape=>shape_rampe
    this%print_t=>print_t_rampe
end subroutine


function shape_rampe(this,time) result(val)
    class (excitation_t),intent(in) :: this
    real(8),intent(in)                  :: time
    real(8)                             :: val(this%dim_mode)

    val=rampe(time,&
              this%dim_mode,&
              this%t_start,&
              this%t_end,&
              this%real_var(1              :  this%dim_mode),&
              this%real_var(this%dim_mode+1:2*this%dim_mode))
end function

function rampe(time,dim_mode,t_start,t_end,val_start,val_end)result(val)
    real(8), intent(in) :: time,t_start,t_end
    integer,intent(in)  :: dim_mode
    real(8), intent(in) :: val_start(dim_mode),val_end(dim_mode)
    real(8)             :: val(dim_mode)

    real(8)             :: t_local

    val=0.0d0
    if ((time>=t_start).and.(time<t_end))then
        t_local=time-t_start
        val=val_start+(val_end-val_start)*t_local/(t_end-t_start)
    endif
end function

subroutine print_t_rampe(this,io)
    class(excitation_t),intent(in)  :: this
    integer,intent(in)              :: io
    character(len=10)   ::  dim_mode

    write(dim_mode,'(I10)') this%dim_mode
    write(io,'(3X,A)') "Time-shape: rampe"
    write(io,'(6X,A)') "Linear shape for I0 to I1 in the time_range"
    write(io,'(6X,A)') "Parameters:"
    write(io,'(9X,A,2(F14.4,A))') "time range:    : [",this%t_start,",",this%t_end," ) fs"
    write(io,'(9X,A,'//dim_mode//'(E16.8))') "start val (I0) : ", this%real_var(1:this%dim_mode)
    write(io,'(9X,A,'//dim_mode//'(E16.8))') "end val   (I1) : ", this%real_var(1+this%dim_mode:this%dim_mode*2)
end subroutine


end module
