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
use m_basic_types, only : vec_point
implicit none
integer, intent(in) :: nim
type(vec_point), intent(in) :: path(:,:)
real(kind=8), intent(out) :: pathlen(nim)
! internal variables
real(kind=8) :: tmp,l1
integer :: i,i_nim
integer :: shape_lattice(2)

pathlen(1) = 0d0
shape_lattice=shape(path)

do i_nim=2,shape_lattice(2)
   tmp = 0d0
   do i=1,shape_lattice(1)

      l1 = calc_ang(path(i,i_nim-1)%w,path(i,i_nim)%w)

      tmp = tmp+l1*l1
   end do
   pathlen(i_nim) = pathlen(i_nim-1) + dsqrt(tmp)
end do

end subroutine the_path





subroutine geodesic_path_one(nim,ni,nf,ax,path,dim_mode)
use m_constants, only : pi
use m_get_random
use m_vector, only : calc_ang
use m_rotation, only : rotation_axis
implicit none
integer, intent(in) :: nim,dim_mode
real(kind=8), intent(in) :: ni(:),nf(:)
real(kind=8), intent(inout) :: ax(:)
real(kind=8), intent(inout) :: path(:,:)
! internal variables
real(kind=8) :: dtheta,theta,angle,tmp,vec(dim_mode),eps=epsilon(angle),pr,norm_local
integer :: i,j

angle = calc_ang(ni,nf)
if (angle<eps) then
   do i = 1,dim_mode
      vec(i) = (nf(i) - ni(i))/(nim-1)
   end do
   path(:,1) = ni(:)
   path(:,nim) = nf(:)
   do i=2,nim-1
      do j=1,dim_mode
         path(j,i) = ni(j) + real((i-1),kind=8)*vec(j)
      end do
      norm_local=norm(path(:,i))
      path(:,i)=path(:,i)/norm_local
   end do

elseif (dabs(angle-pi)<eps) then
      !   write(*,*) 'i am here!'
   pr = 0d0
   do j=1,dim_mode
      pr = pr + ax(j)*ni(j)
   end do
      ax(:) = ax(:) - pr*ni(:)
      tmp = norm(ax)
   do while (tmp<eps)
      pr = 0d0
      do j=1,dim_mode
         ax(j) = 2.0d0*get_rand_classic()
         pr = pr + ax(j)*ni(j)
      end do
         ax(:) = ax(:) - pr*ni(:)
         tmp = norm(ax)
   end do
   norm_local=norm(ax)
   ax=ax/norm_local
   dtheta = pi/real(nim-1)
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
use m_operator_pointer_utils
implicit none
real(kind=8), intent(in) :: amp_rnd
real(kind=8), target, intent(inout) :: path(:,:,:)
! internal variable
integer :: i,i_nim,N_cell,nim,shape_path(3),i_dim,dim_mode,dim_mode_mag,j
real(kind=8) :: norm_local,dumy(3)
real(kind=8), allocatable, dimension(:,:) ::path_one
real(kind=8), allocatable, dimension(:) ::u(:),ax(:),unchanged(:)
real(kind=8), pointer, dimension(:,:,:) :: path_magnetic
logical :: i_magnetic

shape_path=shape(path)
dim_mode=shape_path(1)
N_cell=shape_path(2)
nim=shape_path(3)

call associate_pointer(path_magnetic,path,'magnetic',i_magnetic)
dim_mode_mag=size(path_magnetic(:,1,1))

allocate(path_one(dim_mode_mag,nim),u(dim_mode_mag),ax(dim_mode_mag),unchanged(dim_mode))
path_one=0.0d0
ax = 0.0d0
unchanged=path(:,1,1)

do i_nim = 1,nim
   do i=1,N_cell

   ! part of the modes that will undergo the GNEB
      call geodesic_path_one(nim,path_magnetic(:,i,1),path_magnetic(:,i,nim),ax,path_one,dim_mode_mag)
   ! part of the field and stuff that do not change during optimization
      path(:,i,i_nim)=unchanged
   ! update the part of the path that will change
      do j=1,dim_mode_mag/3
        u(1)=2d0*(get_rand_classic()-1d0)
        u(2)=2d0*(get_rand_classic()-1d0)
        u(3)=2d0*(get_rand_classic()-1d0)
        dumy=path_one( 3*(j-1)+1:3*j ,i_nim)+amp_rnd*u
        norm_local=norm(dumy)
        if (norm_local.gt.1.0d-8) path_magnetic( 3*(j-1)+1:3*j ,i,i_nim)=dumy/norm_local
      enddo

   end do
end do

end subroutine geodesic_path

end module m_path
