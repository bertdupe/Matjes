module m_forces

contains

subroutine forces(tag,field,N_dim,r)
!use m_eval_Beff
use m_basic_types, only : vec_point
use m_derivative, only : calculate_derivative
use m_io_files_utils
use m_convert
use m_vector, only : cross,norm
implicit none

!
!
! Routine that calculates the forces
!
!

integer, intent(in) :: tag,N_dim
real(kind=8), intent(in) :: r(:,:)
type(vec_point), intent(in) :: field(:)
! internals
character(len=50) :: fname
!povray stuff
real(kind=8) :: force(3),dmdr(3,3),B(N_dim),dmdr_int(3,3),volume
integer :: iomp,N_cell,i,io_file

N_cell=size(field)
force=cross(r(1,:),r(2,:))
volume=dot_product(force,r(3,:))

fname=convert('forces_',tag,'.dat')
io_file=open_file_write(fname)

do iomp=1,N_cell

   dmdr=0.0d0
   dmdr_int=0.0d0
   B=0.0d0
!   call calculate_Beff(B,iomp,field,N_dim)
    ERROR STOP "FORCES HAS TO UPDATE CALCULATE_BEFF"
   call calculate_derivative(dmdr_int,iomp)
!
! dmdr_int(:,1) is the derivative along the first unit vector
! dmdr_int(:,2) is the derivative along the second unit vector
! ....
!

   dmdr=matmul(dmdr_int,r)

!
! dmdr(:,1) is the derivative along the cartesian x coordinate
! dmdr(:,2) is the derivative along the cartesian y coordinate
! ....
!

   do i=1,3
      force(i)=dot_product(dmdr(:,i),B(1:3))/volume
   enddo

   write(io_file,'(3(2x,E20.12E3))')(force(i),i=1,3)

enddo

call close_file(fname,io_file)

end subroutine forces

end module m_forces
