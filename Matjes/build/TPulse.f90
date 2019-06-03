module m_TPulse

type Temp_Pulse
    real(kind=8) :: t_cut_pulse,T0,alpha,I_exp,I_over_x,t_exp,t_frac,expo,sigma
    integer :: t_start
end type Temp_Pulse

type(Temp_Pulse) :: TPulse

private
public :: get_parameter_TPulse,update_TPulse
contains

! --------------------------
! initialize the temperature pulse
! --------------------------
subroutine get_parameter_TPulse(io,fname)
use m_io_utils
use m_constants, only : k_b
implicit none
character(len=*), intent(in) :: fname
integer, intent(in) :: io
! internal

! temperature imposee

TPulse%T0=150
TPulse%t_start=25
TPulse%alpha=1.0d0
TPulse%I_exp=137.448
TPulse%I_over_x=286.415
TPulse%t_exp=22.57
TPulse%t_frac=19.62
TPulse%expo=0.47
TPulse%sigma=52.29
TPulse%t_cut_pulse=2500

call get_parameter(io,fname,'T0',TPulse%T0)
call get_parameter(io,fname,'t_start',TPulse%t_start)
call get_parameter(io,fname,'alpha',TPulse%alpha)
call get_parameter(io,fname,'I_exp',TPulse%I_exp)
call get_parameter(io,fname,'I_over_x',TPulse%I_over_x)
call get_parameter(io,fname,'t_exp',TPulse%t_exp)
call get_parameter(io,fname,'t_frac',TPulse%t_frac)
call get_parameter(io,fname,'expo',TPulse%expo)
call get_parameter(io,fname,'sigma',TPulse%sigma)
call get_parameter(io,fname,'t_cut_pulse',TPulse%t_cut_pulse)


TPulse%I_exp=TPulse%I_exp*TPulse%alpha
TPulse%I_over_x=TPulse%I_over_x*TPulse%alpha

end subroutine get_parameter_TPulse


! --------------------------
! update the temperature pulse
! --------------------------
subroutine update_TPulse(time,kt)
use m_constants, only : k_b
implicit none
real(kind=8), intent(in) :: time
real(kind=8), intent(inout) :: kt
! internal
real(kind=8) :: kt1,dumy

kt1=0.0d0

if (time.le.TPulse%t_start) then
    kt1=TPulse%T0
elseif ((time.gt.TPulse%t_start).and.(time.lt.(TPulse%t_start+TPulse%t_cut_pulse))) then
    kt1=TPulse%T0+TPulse%I_exp*exp(-(real(time-TPulse%t_start)-TPulse%t_exp)**2/TPulse%sigma)
elseif (time.ge.(TPulse%t_start+TPulse%t_cut_pulse)) then
    kt1=TPulse%T0+TPulse%I_over_x/(real(time-TPulse%t_start)-TPulse%t_frac)**TPulse%expo
endif

kt=kt1/650.0d0*29.0d0*k_b

end subroutine update_TPulse

end module m_TPulse
