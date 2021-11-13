module m_rotation_matrix
use m_vector, only : norm
use m_constants, only : identity
implicit none
private
public :: rotate_matrix,rotation_matrix

interface rotate_matrix
   module procedure rotate_matrix_real
end interface rotate_matrix

interface rotation_matrix
   module procedure rotation_matrix_real,rotation_matrix_int
end interface rotation_matrix

contains

   subroutine rotate_matrix_real(mat_in,theta,rotation_axis)
   real(8), intent(inout) :: mat_in(:,:)
   real(8), intent(in)    :: theta,rotation_axis(:)

   real(8)                :: u(3),mat(3,3),rot_mat(3,3)
   integer                :: shape_mat_in(2)

   if (abs(theta).lt.1.0d-8) return

   shape_mat_in=shape(mat_in)
   if (size(rotation_axis).ne.3) STOP "ERROR rotation axis should be of size 3"
   if ((shape_mat_in(1).ne.3).or.(shape_mat_in(2).ne.3)) STOP "ERROR the matrix should be 3x3"

   u=rotation_axis/norm(rotation_axis)

   call rotation_matrix(rot_mat,theta,u)
   mat=matmul(mat_in,rot_mat)

   call rotation_matrix(rot_mat,-theta,u)
   mat_in=matmul(rot_mat,mat)

   end subroutine

   subroutine rotation_matrix_real(mat_out,theta,rotation_axis)
   ! subroutine that returns the rotation matrix around the axis rotation_axis
   real(8), intent(out) :: mat_out(3,3)
   real(8), intent(in)  :: theta,rotation_axis(3)

   real(8)              ::  sinang,cosang

   if (abs(theta).lt.1.0d-8) then
      mat_out=identity(1.0d0)
      return
   endif

   sinang = dsin(theta)
   cosang = dcos(theta)

   mat_out(1,1)=rotation_axis(1)**2*(1.0d0-cosang)+cosang
   mat_out(1,2)=rotation_axis(1)*rotation_axis(2)*(1.0d0-cosang)-rotation_axis(3)*sinang
   mat_out(1,3)=rotation_axis(1)*rotation_axis(3)*(1.0d0-cosang)+rotation_axis(2)*sinang

   mat_out(2,1)=rotation_axis(1)*rotation_axis(2)*(1.0d0-cosang)+rotation_axis(3)*sinang
   mat_out(2,2)=rotation_axis(2)**2*(1.0d0-cosang)+cosang
   mat_out(2,3)=rotation_axis(2)*rotation_axis(3)*(1.0d0-cosang)-rotation_axis(1)*sinang

   mat_out(3,1)=rotation_axis(1)*rotation_axis(3)*(1.0d0-cosang)-rotation_axis(2)*sinang
   mat_out(3,2)=rotation_axis(2)*rotation_axis(3)*(1.0d0-cosang)+rotation_axis(1)*sinang
   mat_out(3,3)=rotation_axis(3)**2*(1.0d0-cosang)+cosang
   end subroutine

   subroutine rotation_matrix_int(mat_out,theta,rotation_axis)
   ! subroutine that returns the rotation matrix around the axis rotation_axis
   real(8), intent(out) :: mat_out(3,3)
   real(8), intent(in)  :: theta
   integer, intent(in)  :: rotation_axis(3)

   real(8)              ::  vec(3),vec_tmp(3)

   vec_tmp=real(rotation_axis,8)
   vec=vec_tmp/norm(vec_tmp)

   call rotation_matrix_real(mat_out,theta,vec)

   end subroutine

end module m_rotation_matrix
