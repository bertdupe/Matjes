module m_heaviside
!this contains the heaviside function for the excitation shape
!it is not really a heaviside function, but the name was chosen like that and I won't change it at this point
use m_excitation_shape_base, only: excitation_shape
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit
implicit none
private
public :: update_heavyside
public :: shape_heaviside, heaviside_read_string

contains

subroutine heaviside_read_string(this,str)
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
        write(error_unit,'(A,I3,A)') "The error happens when reading the ", size(this%real_var)," reals after 2 strings, 2 reals, and 1 int"
        STOP
    endif
    this%get_shape=>shape_heaviside
end subroutine

function shape_heaviside(this,time) result(val)
    class (excitation_shape),intent(in) :: this
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
    if ((time.ge.t_start).and.(time.le.t_end)) val=val_start
    if (time.gt.t_end) val=val_end(:)
end function



subroutine update_heavyside(time,field,t_start,t_end,start_value,end_value,counter)
    implicit none
    real(kind=8), intent(in) :: time,t_start,t_end
    real(kind=8), intent(in) :: start_value(:),end_value(:)
    integer, intent(inout) :: counter
    real(kind=8), intent(inout) :: field(:)
    ! internal
    
    if ((time.ge.t_start).and.(time.le.t_end)) field=start_value(:)
    
    
    if (time.gt.t_end) field=end_value(:)

end subroutine update_heavyside

end module 
