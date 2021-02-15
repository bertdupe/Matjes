! module that deals with the velocities of the atoms for the molecular dynamics routine

module m_velocities

contains

subroutine initialize_velocities(V)
use m_io_files_utils
real(8), intent(out) :: V(:,:)
! internal variables
integer :: shape_V(2),io_V,i,j
logical :: file_found=.false.


V=0.0d0

inquire(file='velocity.in', exist=file_found)
if (file_found) then
  shape_V=shape(V)
  io_V=open_file_read('velocity.in')

  do i=1,shape_V(2)
    read(io_V,*) (V(j,i),j=1,shape_V(1))
  enddo

  call close_file('velocity.in',io_V)
else
  write(6,'(a)') "initial speeds are zeros"
endif


end subroutine initialize_velocities


end module m_velocities
