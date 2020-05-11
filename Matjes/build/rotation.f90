module m_rotation

! rotation of vectors
   interface rotate
       module procedure rotate_vector
   end interface rotate

! calculate the rotation axis between 2 vectors
   interface rotation_axis
       module procedure rotation_axis_1D
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

! ==============================================================
!> Calculate axis of rotation based on two 3-vectors
function rotation_axis_1D(m_in,f_in)
use m_vector, only : norm_cross,norm,cross
implicit none
real(kind=8), intent(in) :: m_in(:),f_in(:)
! internal varible
real(kind=8), parameter :: fact=10000000d0
real(kind=8), dimension(3) :: rotation_axis_1D
real(kind=8) :: x(3),tmp,y(3),a,b,c,eps=epsilon(a)


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

end module m_rotation
