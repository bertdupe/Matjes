module m_exc_t_EMwave
use m_exc_t_base, only: excitation_t
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit
implicit none

private
public :: EMwave_read_string

contains


subroutine EMwave_read_string(this,str)
    use m_constants, only: pi
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
    allocate(this%real_var(this%dim_mode*3+2))
    read(str,*,iostat=stat)  dummy_name, shape_name, this%t_start, this%t_end, this%dim_mode, this%real_var
    if(stat/=0)then
        write(error_unit,'(A)') "ERROR, Failed to read excitation shape from line:"
        write(error_unit,'(A)') str
        write(error_unit,'(A,I3,A)') "The error happens when reading the ", size(this%real_var)," reals after 2 strings, 2 reals, and 1 int"
        STOP
    endif

    associate( I0     => this%real_var(1                  :  this%dim_mode*1  ),&
               omega  => this%real_var(1+this%dim_mode*1  :  this%dim_mode*2  ),&
               phi    => this%real_var(1+this%dim_mode*2  :  this%dim_mode*3  ),&
               t0     => this%real_var(1+this%dim_mode*3  :  this%dim_mode*3+1),&
               tau    => this%real_var(1+this%dim_mode*3+1:  this%dim_mode*3+2)) 

        phi=phi*pi              !insert phi in units of pi
        omega=2.0d0*pi/omega    !insert omega as wavelength in 1/fs

        this%t_start=max(this%t_start,t0(1)-10.0d0*tau(1))
        this%t_end  =min(this%t_end  ,t0(1)+10.0d0*tau(1))
    end associate

    this%get_shape=>shape_EMwave
    this%print_t=>print_t_EMwave
end subroutine


function shape_EMwave(this,time) result(val)
    class (excitation_t),intent(in) :: this
    real(8),intent(in)                  :: time
    real(8)                             :: val(this%dim_mode)

    associate( I0     => this%real_var(1                  :  this%dim_mode*1  ),&
               omega  => this%real_var(1+this%dim_mode*1  :  this%dim_mode*2  ),&
               phi    => this%real_var(1+this%dim_mode*2  :  this%dim_mode*3  ),&
               t0     => this%real_var(1+this%dim_mode*3  :  this%dim_mode*3+1),&
               tau    => this%real_var(1+this%dim_mode*3+1:  this%dim_mode*3+2)) 

        val=Emwave(time,&
                   this%dim_mode,&
                   this%t_start,&
                   this%t_end,&
                   I0,omega,phi,t0(1),tau(1))
    end associate
end function

function EMwave(time,dim_mode,t_start,t_end,I0,omega,phi,t0,tau)result(val)
    real(8),intent(in)                      :: time,t_start,t_end
    integer,intent(in)                      :: dim_mode
    real(8),intent(in),dimension(dim_mode)  :: I0,omega,phi
    real(8),intent(in)                      :: t0, tau
    real(8)                                 :: val(dim_mode)
    real(8)     :: tmp

    val=0.d0
    if ((time>=t_start).and.(time<t_end))then
        tmp=time-t0
        tmp=tmp/tau
        tmp=-tmp*tmp
        if(tmp>-200.0d0)then    !prevent underflow
            tmp=exp(tmp)
            val=I0 *cos(omega*time+phi)*tmp
        endif
    endif
end function

subroutine print_t_EMwave(this,io)
    use m_constants, only: pi
    class(excitation_t),intent(in)  :: this
    integer,intent(in)              :: io
    character(len=10)   ::  dim_mode

    write(dim_mode,'(I10)') this%dim_mode
    associate( I0     => this%real_var(1                  :  this%dim_mode*1  ),&
               omega  => this%real_var(1+this%dim_mode*1  :  this%dim_mode*2  ),&
               phi    => this%real_var(1+this%dim_mode*2  :  this%dim_mode*3  ),&
               t0     => this%real_var(1+this%dim_mode*3  :  this%dim_mode*3+1),&
               tau    => this%real_var(1+this%dim_mode*3+1:  this%dim_mode*3+2)) 
    write(io,'(3X,A)') "Time-shape: EMwave"
    write(io,'(6X,A)') "Wave with mulitplicative gaussian shape"
    write(io,'(6X,A)') "Equation: I0*cos(w*t+phi)*exp(-(t-t0)^2/tau^2)"
    write(io,'(6X,A)') "Parameters:"
    write(io,'(9X,A,2(F14.4,A))') "time range:    : [",this%t_start,",",this%t_end," ) fs"
    write(io,'(9X,A,'//dim_mode//'(E16.8))') "amplitude (I0) : ", I0
    write(io,'(9X,A,'//dim_mode//'(F16.8),A)') "frequency  (w) : ", omega/2.0d0/pi, " in 2pi/fs"
    write(io,'(9X,A,'//dim_mode//'(F16.4),A)') "Phase    (phi) : ", phi/pi , " in pi"
    write(io,'(9X,A,(F16.4),A)') "center    (t0) : ", t0(1) , " in fs"
    write(io,'(9X,A,(F16.4),A)') "decay    (tau) : ", tau(1) , " in fs"
    end associate
end subroutine

end module
