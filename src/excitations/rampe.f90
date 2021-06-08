module m_rampe
use m_convert
use m_excitation_shape_base, only: excitation_shape
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit

private
public :: update_rampe
public :: shape_rampe, rampe_read_string

contains


subroutine rampe_read_string(this,str)
    class(excitation_shape),intent(inout)   :: this
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
end subroutine


function shape_rampe(this,time) result(val)
    class (excitation_shape),intent(in) :: this
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
    if ((time.ge.t_start).and.(time.le.t_end))then
        t_local=time-t_start
        val=(val_end-val_start)*t_local/(t_end-t_start)
    endif
end function


subroutine update_rampe(time,field,t_start,t_end,start_value,end_value,counter)
implicit none
real(kind=8), intent(in) :: time,t_start,t_end
real(kind=8), intent(in) :: end_value(:),start_value(:)
integer, intent(inout) :: counter
real(kind=8), intent(inout) :: field(:)
! internal
character(len=30) :: form
integer :: n

n=size(field)

if ((time.ge.t_start).and.(time.le.t_end)) then

   field=field+(end_value-start_value)/real(t_end-t_start,8)

   form=convert('(a,',n,'f14.6)')
   write(6,form) 'field value ',field
endif

end subroutine update_rampe

end module m_rampe
