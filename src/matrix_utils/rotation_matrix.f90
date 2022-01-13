module m_rotation_matrix
use m_vector, only : norm
use m_constants, only : identity
implicit none
private
public :: rotate_matrix,check_rotate_matrix

interface rotate_matrix
   module procedure rotation_matrix_int,rotation_matrix_real_sym_op,rotate_matrix_real
end interface rotate_matrix

contains

   subroutine check_rotate_matrix(theta,rotation_axis,bound_input,vec)
   real(8), intent(inout) :: theta
   real(8), intent(in)    :: rotation_axis(3),bound_input(3),vec(3)

   real(8) :: u(3),rot_mat(3,3),vec_tmp(3)

   u=rotation_axis/norm(rotation_axis)

   call rotation_matrix_real(rot_mat,theta,u)

   vec_tmp=matmul(rot_mat,bound_input)

   if (norm(vec_tmp-vec).gt.1.0d-5) STOP 'angle not correct in check_rotate_matrix'

   end subroutine

   subroutine rotation_matrix_real_sym_op(mat_in,sym_op)
   real(8), intent(inout) :: mat_in(:,:)
   real(8), intent(in)    :: sym_op(:,:)

   real(8)                :: mat(3,3)
   integer                :: shape_mat_in(2)

   shape_mat_in=shape(mat_in)
   if ((shape_mat_in(1).ne.3).or.(shape_mat_in(2).ne.3)) STOP "ERROR the matrix should be 3x3"

   mat=matmul(mat_in,sym_op)
   mat_in=matmul(transpose(sym_op),mat)

   end subroutine

   subroutine rotate_matrix_real(mat_out,mat_in,theta,rotation_axis)
   real(8), intent(out)   :: mat_out(:,:)
   real(8), intent(in)    :: mat_in(:,:)
   real(8), intent(in)    :: theta,rotation_axis(:)

   real(8)                :: u(3),mat(3,3),rot_mat(3,3)
   integer                :: shape_mat_in(2)

   if (abs(theta).lt.1.0d-8) then
       mat_out=mat_in
       return
   endif

   shape_mat_in=shape(mat_in)
   if (size(rotation_axis).ne.3) STOP "ERROR rotation axis should be of size 3"
   if ((shape_mat_in(1).ne.3).or.(shape_mat_in(2).ne.3)) STOP "ERROR the matrix should be 3x3"

   u=rotation_axis/norm(rotation_axis)

   call rotation_matrix_real(rot_mat,theta,u)
   mat=matmul(mat_in,rot_mat)
   mat_out=matmul(transpose(rot_mat),mat)

   end subroutine

   subroutine rotation_matrix_real(mat_out,theta,rotation_axis)
   ! subroutine that returns the rotation matrix around the axis rotation_axis
   real(8), intent(out) :: mat_out(3,3)
   real(8), intent(in)  :: theta,rotation_axis(3)

   real(8)              ::  sinang,cosang

   if ((norm(rotation_axis)-1.0d0).gt.1.0d-8) STOP 'rotation axis should be normalized to 1 in rotation_matrix_real'

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
