module m_update_time
!!!!!!!!!!!!!!!!!!!!!
! module that update the time as a function of the effective magnetic field
!!!!!!!!!!!!!!!!!!!!!
logical :: i_update_time
real(kind=8) :: discretization
interface update_time
  module procedure update_time_only_B
    
end interface

private
public :: update_time,init_update_time

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

subroutine update_time_only_B(timestep,B,damping)
use m_constants, only : hbar,pi
use m_vector, only : norm
implicit none
real(kind=8), intent(inout) :: timestep
real(kind=8), intent(in) :: damping
real(kind=8),dimension(:,:),intent(in) :: B
! internal
integer :: size_B
real(kind=8) :: max_B,dumy_B,timestep_backup
integer :: i

if (.not.i_update_time) return
max_B=00.0d0
size_B=size(B,2)
timestep_backup=timestep

do i=1,size_B
   dumy_B=norm(B(:,i))

!!!!!!!!!!!!!!!!!!!!!!!!
! this leads to too large time steps in case of laser pulses
!   if ((dumy_B+dumy_BT).gt.max_B) max_B=dumy_B/(1.0d0+damping**2)+dumy_BT*sqrt(damping/(1.0d0+damping**2))
!!!!!!!!!!!!!!!!!!!!!!!!
   if ((dumy_B).gt.max_B) max_B=dumy_B

enddo

if (max_B.gt.1.0d-8) then
   timestep=2.0d0*pi*hbar/max_B/discretization
else
   stop 'error in update_time'
endif

end subroutine 


end module m_update_time
