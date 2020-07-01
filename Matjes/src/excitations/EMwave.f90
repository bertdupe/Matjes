module m_EMwave
!!!!!!!!!!!!!!!!!!!!!!!!!
! module that deals with the electromagnetic excitations
! with spatial and temporal dependences
!
! I in W/cm^2 as input
! E in V/nm in the code
!
!!!!!!!!!!!!!!!!!!!!!!!!!

type EMwave
    real(kind=8) :: omega_l,Tau,E_0,t_0,t_start,t_cut
    !integer :: t_start,t_cut
end type EMwave

type(EMwave) :: EM_Pulse

private
public :: get_parameter_EMwave,update_EMwave

contains

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
!EM_Pulse%omega_l=pi(2.0d0)*c/lambda_l
EM_Pulse%omega_l=pi(2.0d0)/lambda_l

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
EM_Pulse%E_0=I_0*alpha

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
	!if (field(3).lt.1.0d-8) field=0.0d0
else !before and after wave
	field(3)=0.0d0
endif


end subroutine update_EMwave

end module m_EMwave
