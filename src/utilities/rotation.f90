module m_rotation
implicit none

! rotation of vectors
   interface rotate
       module procedure rotate_vector
       module procedure rotate_vector_arr
   end interface rotate

! calculate the rotation axis between 2 vectors
   interface rotation_axis
       module procedure rotation_axis_1D
       module procedure rotation_axis_arr
   end interface rotation_axis

private
public :: rotate,rotation_axis
contains

! ==============================================================
!> Rotate 3-vector v_in by angle ang around axis ax
subroutine rotate_vector(v_in,ax,ang,v_out)
implicit none
real(kind=8), intent(in) :: v_in(3),ax(3),ang
real(kind=8), intent(out) :: v_out(3)
! internal variables
real(kind=8) :: tmp,sinang,cosang

tmp = dot_product(ax,v_in)

sinang = dsin(ang)
cosang = dcos(ang)

v_out(1) = v_in(1)*cosang + sinang*(ax(2)*v_in(3)-ax(3)*v_in(2))+(1d0-cosang)*tmp*ax(1)
v_out(2) = v_in(2)*cosang - sinang*(ax(1)*v_in(3)-ax(3)*v_in(1))+(1d0-cosang)*tmp*ax(2)
v_out(3) = v_in(3)*cosang + sinang*(ax(1)*v_in(2)-ax(2)*v_in(1))+(1d0-cosang)*tmp*ax(3)

end subroutine rotate_vector

!> Rotate 3-vector v_in by angle ang around axis ax
pure subroutine rotate_vector_arr(v_in,ax,ang,v_out)
use  m_vector, only: cross
!!$ use omp_lib
implicit none
real(8), intent(in)         :: v_in(:,:),ax(:,:),ang(:)
real(8), intent(out)        :: v_out(size(v_in,1),size(v_in,2))
! internal variables
real(8)     :: tmp(size(v_in,2)),sinang(size(v_in,2)),cosang(size(v_in,2))
integer     :: i

tmp = sum(ax*v_in,1)
sinang = dsin(ang)
cosang = dcos(ang)
do i=1,size(v_in,2)
    v_out(1,i)=v_in(1,i)*cosang(i)+sinang(i)*(ax(3,i)*v_in(2,i)-ax(2,i)*v_in(3,i))+(1.0d0-cosang(i))*tmp(i)*ax(1,i)
    v_out(2,i)=v_in(2,i)*cosang(i)+sinang(i)*(ax(1,i)*v_in(3,i)-ax(3,i)*v_in(1,i))+(1.0d0-cosang(i))*tmp(i)*ax(2,i)
    v_out(3,i)=v_in(3,i)*cosang(i)+sinang(i)*(ax(2,i)*v_in(1,i)-ax(1,i)*v_in(2,i))+(1.0d0-cosang(i))*tmp(i)*ax(3,i)
end do
end subroutine 


! ==============================================================
!> Calculate axis of rotation based on two 3-vectors
function rotation_axis_1D(m_in,f_in)
use m_vector, only : norm_cross,norm,cross
implicit none
real(8), intent(in) :: m_in(:),f_in(:)
! internal varible
real(8), parameter :: fact=10000000d0
real(8), dimension(3) :: rotation_axis_1D
real(8) :: x(3),tmp,y(3),a,b,c,eps=epsilon(a)


tmp = norm(f_in)
if (tmp<eps) then
    a = f_in(1)*fact
    b = f_in(2)*fact
    c = f_in(3)*fact
else
    a = f_in(1)
    b = f_in(2)
    c = f_in(3)
end if


y(1) = dsign(a,f_in(1))
y(2) = dsign(b,f_in(2))
y(3) = dsign(c,f_in(3))

x=cross(m_in,y)
tmp = norm_cross(m_in,y,1,3)
if (tmp.lt.eps) then
  rotation_axis_1D = 0.0d0
else
  rotation_axis_1D = x/tmp
endif

end function rotation_axis_1D

! ==============================================================
!> Calculate axis of rotation based on two 3-vectors in array
function rotation_axis_arr(m_in,f_in)result(res)
    use m_vector, only : norm_cross,norm,cross,normalize,vec_normalize
    real(8),intent(in),target,contiguous   :: m_in(:,:),f_in(:,:)
    real(8),target                         :: res(3,size(m_in,2))
    real(8),pointer         :: m_flat(:),res_flat(:),tmp_flat(:)
    ! internal varible
    real(8),target          :: tmp(3,size(m_in,2))
    real(8)                 :: eps=epsilon(f_in)
  
    m_flat(1:size(m_in))=>m_in
    res_flat(1:size(res))=>res
    tmp_flat(1:size(tmp))=>tmp

    tmp=f_in
    Call normalize(tmp,1.0d-30)
    res_flat=cross(m_flat,tmp_flat,1,size(m_in))
    Call normalize(res,eps)
end function 


end module m_rotation
