module m_torques
use m_derived_types, only : vec_point

type torque
     logical :: on
     real(kind=8) :: torque_AFL,torque_FL
     real(kind=8) :: direction(3)
     type(vec_point), allocatable :: distribution(:)
end type torque

type(torque) :: Slonczewski,SOT

private
public :: get_torques, update_DS, update_B

contains

!!!!!!!!!!!!!!!!!
!
! initialize the different torque to add them up to the neighbouring field
!
!!!!!!!!!!!!!!!!!
subroutine get_torques(fname)
use m_io_utils
use m_io_files_utils
implicit none
character(len=*), intent(in) :: fname
! internal
integer :: io,N_site,i
real(kind=8), allocatable :: pos(:,:)

! temperature imposee
Slonczewski%direction=0.0d0
Slonczewski%on=.false.
Slonczewski%torque_AFL=0.0d0
Slonczewski%torque_FL=0.0d0

SOT%direction=0.0d0
SOT%on=.false.
SOT%torque_AFL=0.0d0
SOT%torque_FL=0.0d0

io=open_file_read(fname)

call get_parameter(io,fname,'Slonczewski_torque',Slonczewski%on)
if (Slonczewski%on) then
  write(6,'(a)') 'The Slonczewski torque is on'
  call get_parameter(io,fname,'Slonczewski_direction',3,Slonczewski%direction)

! get the form of the torques from the positions that are printed in the file
  N_site=get_lines('positions.dat')
  allocate(pos(3,N_site))

  allocate(Slonczewski%distribution(N_site))

  do i=1,N_site
    nullify(Slonczewski%distribution(i)%w)
  enddo

  deallocate(pos)

  call get_parameter(io,fname,'Slonczewski_AFL',Slonczewski%torque_AFL)
  call get_parameter(io,fname,'Slonczewski_FL',Slonczewski%torque_FL)
endif

call get_parameter(io,fname,'SOT_torque',SOT%on)
if (SOT%on) then
  write(6,'(a)') 'The SOT torque is on'
  call get_parameter(io,fname,'SOT_direction',3,SOT%direction,1.0d0)

  call get_parameter(io,fname,'SOT_AFL',SOT%torque_AFL)
  call get_parameter(io,fname,'SOT_FL',SOT%torque_FL)
endif

call close_file(fname,io)

end subroutine

!!!!!!!!!!!!!!!!!
!
! update the different effective field B_eff with the torques contribution
!
!!!!!!!!!!!!!!!!!
subroutine update_DS(S,damping,ds)
use m_vector, only : cross
implicit none
real(kind=8), intent(in) :: damping,S(:)
real(kind=8), intent(inout) :: ds(:)
! internal
real(kind=8) :: steptor(3)

steptor=0.0d0

!if (Slonczewski%on) then
!     steptor=cross(S,Slonczewski%distribution(i)%w)
!     ds=ds+Slonczewski%torque_FL*(1.0d0-damping*Slonczewski%torque_AFL)*steptor
!     ds=ds+Slonczewski%torque_FL*(Slonczewski%torque_AFL+damping)*cross(S,steptor)
!endif

if (SOT%on) then
   steptor=cross(S,SOT%direction)
   ds=ds+SOT%torque_FL*(1.0d0-damping*SOT%torque_AFL)*steptor
   ds=ds+SOT%torque_FL*(SOT%torque_AFL+damping)*cross(S,steptor)
endif

end subroutine

!!!!!!!!!!!!!!!!!
!
! prepare B to write everything under the Y=cross(X,B)
!
!!!!!!!!!!!!!!!!!
subroutine update_B(S,damping,B)
use m_vector, only : cross
implicit none
real(kind=8), intent(in) :: damping,S(:)
real(kind=8), intent(inout) :: B(:)
! internal
real(kind=8) :: steptor(3)

steptor=0.0d0

!if (Slonczewski%on) then
!     steptor=cross(S,Slonczewski%distribution(i)%w)
!     B=B+Slonczewski%torque_FL*(1.0d0-damping*Slonczewski%torque_AFL)*Slonczewski%distribution(i)%w
!     B=B+Slonczewski%torque_FL*(Slonczewski%torque_AFL+damping)*steptor
!endif

if (SOT%on) then
   steptor=cross(S,SOT%direction)
   B=B+SOT%torque_FL*(1.0d0-damping*SOT%torque_AFL)*SOT%direction
   B=B+SOT%torque_FL*(SOT%torque_AFL+damping)*steptor
endif

end subroutine
end module m_torques
