module m_path
!
! files that contains the routine to calculate the path for the GNEB
!
!
!
use m_vector, only : norm

private
public :: the_path,geodesic_path

contains

!> Calculate poly-geodesic length of the path
subroutine the_path(nim,path,pathlen)
use m_vector, only : calc_ang
implicit none
integer, intent(in) :: nim
real(kind=8), intent(in) :: path(:,:,:)
real(kind=8), intent(out) :: pathlen(nim)
! internal variables
real(kind=8) :: tmp,l1
integer :: i,i_nim
integer :: shape_lattice(3)

pathlen(1) = 0d0
shape_lattice=shape(path)

do i_nim=2,shape_lattice(3)
   tmp = 0d0
   do i=1,shape_lattice(2)

      l1 = calc_ang(path(:,i,i_nim-1),path(:,i,i_nim))

      tmp = tmp+l1*l1
   end do
   pathlen(i_nim) = pathlen(i_nim-1) + dsqrt(tmp)
end do

end subroutine the_path





subroutine geodesic_path_one(nim,ni,nf,ax,path)
use m_constants, only : pi
use m_get_random
use m_vector, only : calc_ang
use m_rotation, only : rotation_axis
implicit none
integer, intent(in) :: nim
real(kind=8), intent(in) :: ni(3),nf(3)
real(kind=8), intent(inout) :: ax(3)
real(kind=8), intent(inout) :: path(:,:)
! internal variables
real(kind=8) :: dtheta,theta,angle,tmp,vec(3),eps=epsilon(angle),pr,norm_local
integer :: i,j
type(mtprng_state) :: state

angle = calc_ang(ni,nf)
if (angle<eps) then
   do i = 1,3
      vec(i) = (nf(i) - ni(i))/(nim-1)
   end do
   path(:,1) = ni(:)
   path(:,nim) = nf(:)
   do i=2,nim-1
      do j=1,3
         path(j,i) = ni(j) + real((i-1),kind=8)*vec(j)
      end do
      norm_local=norm(path(:,i))
      path(:,i)=path(:,i)/norm_local
   end do

elseif (dabs(angle-pi(1.0d0))<eps) then
      !   write(*,*) 'i am here!'
   pr = 0d0
   do j=1,3
      pr = pr + ax(j)*ni(j)
   end do
      ax(:) = ax(:) - pr*ni(:)
      tmp = norm(ax)
   do while (tmp<eps)
      pr = 0d0
      do j=1,3
         ax(j) = 2.0d0*get_rand_classic(state)
         pr = pr + ax(j)*ni(j)
      end do
         ax(:) = ax(:) - pr*ni(:)
         tmp = norm(ax)
   end do
   norm_local=norm(ax)
   ax=ax/norm_local
   dtheta = pi(1.0d0)/real(nim-1)
   path(:,1) = ni(:)
   path(:,nim) = nf(:)

   do i=2,nim-1
      theta = (i-1)*dtheta
      path(1,i) = ni(1)*dcos(theta) + dsin(theta)*(ax(2)*ni(3)-ax(3)*ni(2))
      path(2,i) = ni(2)*dcos(theta) - dsin(theta)*(ax(1)*ni(3)-ax(3)*ni(1))
      path(3,i) = ni(3)*dcos(theta) + dsin(theta)*(ax(1)*ni(2)-ax(2)*ni(1))

      norm_local=norm(path(:,i))
      path(:,i)=path(:,i)/norm_local

   end do

else
   ax = rotation_axis(ni,nf)

   dtheta = angle/(nim-1)

   path(:,1) = ni(:)
   path(:,nim) = nf(:)



   do i=2,nim-1
      theta = (i-1)*dtheta
      path(1,i) = ni(1)*dcos(theta) + dsin(theta)*(ax(2)*ni(3)-ax(3)*ni(2))
      path(2,i) = ni(2)*dcos(theta) - dsin(theta)*(ax(1)*ni(3)-ax(3)*ni(1))
      path(3,i) = ni(3)*dcos(theta) + dsin(theta)*(ax(1)*ni(2)-ax(2)*ni(1))

      norm_local=norm(path(:,i))
      path(:,i)=path(:,i)/norm_local

   end do

end if

end subroutine geodesic_path_one


!
! Calculate the geodesic path between all the configurations
!
! I have removed the first configurations spini and spinf but we can put it back

subroutine geodesic_path(amp_rnd,path)
use m_get_random
use m_basic_types, only : vec
implicit none
real(kind=8), intent(in) :: amp_rnd
real(kind=8), intent(inout) :: path(:,:,:)
! internal variable
integer :: i,i_nim,N_cell,nim,shape_path(3)
real(kind=8) :: u(3),ax(3),norm_local
real(kind=8), allocatable, dimension(:,:) ::path_one
type(mtprng_state) :: state

ax = 0d0
shape_path=shape(path)
N_cell=shape_path(2)
nim=shape_path(3)

allocate(path_one(3,nim))
path_one=0.0d0

do i_nim = 1,nim
   do i=1,N_cell

      call geodesic_path_one(nim,path(:,i,1),path(:,i,nim),ax,path_one)

      path(:,i,i_nim) = path_one(:,i_nim)

      u(1)=2d0*(get_rand_classic(state)-1d0)
      u(2)=2d0*(get_rand_classic(state)-1d0)
      u(3)=2d0*(get_rand_classic(state)-1d0)
      path(:,i,i_nim) = path(:,i,i_nim)+amp_rnd*u
      norm_local=norm(path(:,i,i_nim))
      path(:,i,i_nim)=path(:,i,i_nim)/norm_local

   end do
end do

end subroutine geodesic_path

end module m_path
