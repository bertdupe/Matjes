module m_exc_t_gauss
!this contains the gauss function for the excitation shape
use m_exc_t_base, only: excitation_t
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit
implicit none
private
public :: gauss_read_string

contains

subroutine gauss_read_string(this,str)
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
        write(error_unit,'(AI3,A)') "The error happens when reading the ", size(this%real_var)," reals after 2 strings,2 reals and 1 int"
        STOP
    endif
    this%get_shape=>shape_gauss

    associate( I0     => this%real_var(1                :  this%dim_mode  ),&
               t0     => this%real_var(1+this%dim_mode  :  this%dim_mode+1),&
               tau    => this%real_var(1+this%dim_mode+1:  this%dim_mode+2)) 

        this%t_start=max(this%t_start,t0(1)-10.0d0*tau(1))
        this%t_end  =min(this%t_end  ,t0(1)+10.0d0*tau(1))
    end associate
end subroutine

function shape_gauss(this,time) result(val)
    class (excitation_t),intent(in)     :: this
    real(8),intent(in)                  :: time
    real(8)                             :: val(this%dim_mode)

    associate( I0     => this%real_var(1                  :  this%dim_mode*1  ),&
               t0     => this%real_var(1+this%dim_mode*2  :  this%dim_mode*2+1),&
               tau    => this%real_var(1+this%dim_mode*2+1:  this%dim_mode*2+2)) 
    val=gauss(time,&
                  this%dim_mode,&
                  this%t_start,&
                  this%t_end,&
                  I0,t0(1),tau(1))
    end associate
end function

function gauss(time,dim_mode,t_start,t_end,I0,t0,tau)result(val)
    real(8),intent(in)  :: time,t_start,t_end
    integer,intent(in)  :: dim_mode
    real(8),intent(in)  :: I0(dim_mode)
    real(8),intent(in)  :: t0,tau
    real(8)             :: val(dim_mode)

    real(8)     :: tmp

    val=0.0d0
    if ((time>=t_start).and.(time<t_end))then
        tmp=time-t0
        tmp=tmp/tau
        tmp=-0.5d0*tmp*tmp
        tmp=exp(tmp)
        val=val+I0*tmp
    endif
end function
end module 
