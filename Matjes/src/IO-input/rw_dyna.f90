subroutine rw_dyna(N,my_lattice,timestep,damping,Efreq,duration)
use m_constants
use m_derived_types
use m_io_files_utils
use m_io_utils
implicit none
type(lattice), intent(in) :: my_lattice
integer, intent(in) :: N(3)
real(kind=8), intent(out) :: timestep,damping
integer, intent(out) :: duration,Efreq
! internal
integer :: io
real(kind=8) :: net(3,3)
logical :: stmtemp,Ffield,i_Efield

net=my_lattice%areal
Efreq=1
timestep=1.0d0
damping=0.0d0

io=open_file_read('input')

call get_parameter(io,'input','timestep',timestep)
call get_parameter(io,'input','Efreq',Efreq)
call get_parameter(io,'input','duration',duration)
call get_parameter(io,'input','STMtemp',stmtemp)
call get_parameter(io,'input','damping',damping)

call get_parameter(io,'input','Ffield',Ffield)
call get_parameter(io,'input','Efield',i_Efield)

call close_file('input',io)

!if (stmtorque) then
!
!   write(*,*) 'SPSTM torque set up'
!   xmed=dble(N(1))/2.0d0*net(1,1)+dble(N(2))/2.0d0*net(2,1)
!   ymed=dble(N(1))/2.0d0*net(1,2)+dble(N(2))/2.0d0*net(2,2)
!
!   do i_x=1,N(1)
!     do i_y=1,N(2)
!       do i_z=1,N(3)
!
!       htor(i_x,i_y,i_z)=dexp(-dsqrt((xmed-my_lattice%l_modes(i_x,i_y,i_z,1)%w(1))**2+ &
!           (ymed-my_lattice%l_modes(i_x,i_y,i_z,1)%w(2))**2+h**2)/ri)
!
!       enddo
!     enddo
!   enddo
!endif
end subroutine rw_dyna
