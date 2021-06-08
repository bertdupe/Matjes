module m_EMwave
!!!!!!!!!!!!!!!!!!!!!!!!!
! module that deals with the electromagnetic excitations
! with spatial and temporal dependences
!
! I in W/cm^2 as input
! E in V/nm in the code
!
!!!!!!!!!!!!!!!!!!!!!!!!!
use m_excitation_shape_base, only: excitation_shape
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit
implicit none

type EMwave_t
    real(kind=8) :: omega_l,Tau,E_0,t_0,t_start,t_cut
    !integer :: t_start,t_cut
end type 

type(EMwave_t) :: EM_Pulse

private
public :: get_parameter_EMwave,update_EMwave
public :: EMwave_read_string

contains


subroutine EMwave_read_string(this,str)
    use m_constants, only: pi
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
    allocate(this%real_var(this%dim_mode*3+2))
    read(str,*,iostat=stat)  dummy_name, shape_name, this%t_start, this%t_end, this%dim_mode, this%real_var
    if(stat/=0)then
        write(error_unit,'(A)') "ERROR, Failed to read excitation shape from line:"
        write(error_unit,'(A)') str
        write(error_unit,'(AI3,A)') "The error happens when reading the ", size(this%real_var)," reals after 2 strings, 2 reals, and 1 int"
        STOP
    endif

    associate( I0     => this%real_var(1                  :  this%dim_mode*1  ),&
               omega  => this%real_var(1+this%dim_mode*1  :  this%dim_mode*2  ),&
               phi    => this%real_var(1+this%dim_mode*2  :  this%dim_mode*3  ),&
               t0     => this%real_var(1+this%dim_mode*3  :  this%dim_mode*3+1),&
               tau    => this%real_var(1+this%dim_mode*3+1:  this%dim_mode*3+2)) 

        phi=phi*pi              !insert phi in units of pi
        omega=2.0d0*pi/omega    !insert omega as wavelength in 1/fs

    end associate

    this%get_shape=>shape_EMwave
end subroutine


function shape_EMwave(this,time) result(val)
    class (excitation_shape),intent(in) :: this
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
    if ((time.ge.t_start).and.(time.le.t_end))then
        tmp=time-t0
        tmp=tmp/tau
        tmp=-tmp*tmp
        if(tmp>-200.0d0)then    !prevent underflow
            tmp=exp(tmp)
            val=I0 *cos(omega*time+phi)*tmp
        endif
    endif
end function




!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize parameters
!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_parameter_EMwave(io,fname)
use m_io_utils
use m_constants, only : epsilon_0,pi,c
implicit none
character(len=*), intent(in) :: fname
integer, intent(in) :: io
! internal
! conversion factor from the intensity to the field
real(kind=8), parameter :: alpha=2.744923727d-6
real(kind=8) :: I_0,lambda_l

EM_Pulse%t_start=-200.0d0
EM_Pulse%t_cut=-800.0d0
lambda_l=25.0d0
EM_Pulse%Tau=100.0d0
I_0=1.0d7
EM_Pulse%t_0=500.0d0

call get_parameter(io,fname,'t_start',EM_Pulse%t_start)
call get_parameter(io,fname,'t_cut',EM_Pulse%t_cut)
call get_parameter(io,fname,'lambda_l',lambda_l)
call get_parameter(io,fname,'Tau',EM_Pulse%Tau)
call get_parameter(io,fname,'I_0',I_0)
call get_parameter(io,fname,'t_0',EM_Pulse%t_0)

! en fs-1
!EM_Pulse%omega_l=2.0d0*pi*c/lambda_l
EM_Pulse%omega_l=2.0d0*pi/lambda_l

if (EM_Pulse%t_start.lt.0) then
   write(6,'(a)') 't_start is negative or not read in input'
   EM_Pulse%t_start=-EM_Pulse%t_start
else
   EM_Pulse%t_start=int(EM_Pulse%t_0-3*EM_Pulse%Tau)
endif

if (EM_Pulse%t_cut.lt.0) then
   write(6,'(a)') 't_cut is negative or not read in input'
   EM_Pulse%t_cut=-EM_Pulse%t_cut
else
   EM_Pulse%t_cut=int(EM_Pulse%t_0+3*EM_Pulse%Tau)
endif

! the electric field is now in V/nm
EM_Pulse%E_0=I_0!*alpha
!no alpha trying to do a Bfield wave

end subroutine get_parameter_EMwave

!!!!!!!!!!!!!!!!!!!!!!!!!
! update EMwave
!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine update_EMwave(time,field)
use m_io_utils
use m_constants, only : epsilon_0
implicit none
real(kind=8), intent(in) :: time
real(kind=8), intent(inout) :: field(:)
! internal

field=0.0d0

if ( (time.ge.EM_Pulse%t_start).and.(time.le.EM_Pulse%t_cut) ) then !during wave
	field(3)=EM_Pulse%E_0 *cos(EM_Pulse%omega_l*time) * exp(-((time-EM_Pulse%t_0)/EM_Pulse%Tau)**2)
else !before and after wave
	field(3)=0.0d0
endif


end subroutine update_EMwave

end module m_EMwave
