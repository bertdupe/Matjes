module m_update_time
!!!!!!!!!!!!!!!!!!!!!
! module that update the time as a function of the effective magnetic field
!!!!!!!!!!!!!!!!!!!!!
logical :: i_update_time
real(kind=8) :: discretization

private
public :: update_time,init_update_time,get_dt_LLG,multiply_2,multiply_1

contains

!!!!!!!!!!!!!!!!!!!!!
! initialize the time update
!!!!!!!!!!!!!!!!!!!!!
subroutine init_update_time(fname)
use m_io_utils
use m_io_files_utils
implicit none
character(len=*), intent(in) :: fname
! internal
integer :: io

i_update_time=.false.
discretization=15.0d0

io=open_file_read(fname)

call get_parameter(io,fname,'update_time',i_update_time)
if (.not.i_update_time) return
call get_parameter(io,fname,'discretization',discretization)

call close_file(fname,io)

end subroutine init_update_time

!!!!!!!!!!!!!!!!!!!!!
! update the time depending on the effective field
!!!!!!!!!!!!!!!!!!!!!
subroutine update_time(timestep,B,BT,damping)
use m_constants, only : hbar,pi
use m_vector, only : norm
implicit none
real(kind=8), intent(inout) :: timestep
real(kind=8), intent(in) :: damping
real(kind=8),dimension(:,:),intent(in) :: B,BT
! internal
integer :: size_B,size_BT
real(kind=8) :: max_B,dumy_B,dumy_BT,timestep_backup
integer :: i

if (.not.i_update_time) return
max_B=00.0d0
dumy_BT=0.0d0
size_B=size(B,2)
size_BT=size(BT,2)
timestep_backup=timestep

if (size_B.ne.size_BT) stop 'error in update_time'

do i=1,size_B
   dumy_B=norm(B(:,i))
   dumy_BT=dumy_BT+norm(BT(:,i))

!!!!!!!!!!!!!!!!!!!!!!!!
! this leads to too large time steps in case of laser pulses
!   if ((dumy_B+dumy_BT).gt.max_B) max_B=dumy_B/(1.0d0+damping**2)+dumy_BT*sqrt(damping/(1.0d0+damping**2))
!!!!!!!!!!!!!!!!!!!!!!!!
   if ((dumy_B).gt.max_B) max_B=dumy_B

enddo

if (max_B.gt.1.0d-8) then
   dumy_BT=dumy_BT/real(size_BT)
   timestep=pi(2.0d0)*hbar/(max_B/(1.0d0+damping**2)+dumy_BT*sqrt(damping/(1.0d0+damping**2)))/discretization
else
   stop 'error in update_time'
endif

!write(6,'(a,f8.4)') 'The old timestep is', timestep_backup
!write(6,'(a,f8.4)') 'The new timestep is', timestep

end subroutine update_time

real(kind=8) function get_dt_LLG(timestep,damping)
use m_constants, only : hbar
implicit none
real(kind=8) :: timestep,damping

get_dt_LLG=timestep/hbar/(1.0d0+damping**2)

end function

real(kind=8) function multiply_1(dt)
implicit none
real(kind=8), intent(in) :: dt

multiply_1=1.0d0*dt

end function

real(kind=8) function multiply_2(dt)
implicit none
real(kind=8), intent(in) :: dt

multiply_2=2.0d0*dt

end function

end module m_update_time
