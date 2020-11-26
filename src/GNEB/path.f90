module m_path
!
! files that contains the routine to calculate the path for the GNEB
!
!
!
use m_type_lattice,only: lattice
use m_vector, only : norm
implicit none

private
public :: the_path,geodesic_path

contains

!> Calculate poly-geodesic length of the path
subroutine the_path(images,pathlen)
    use m_vector, only : calc_ang
    implicit none
    type(lattice), intent(in)   :: images(:)
    real(8), intent(out)        :: pathlen(size(images))
    ! internal variables
    real(8) :: tmp,l1
    integer :: i,i_nim
    
    pathlen(1) = 0d0
    
    do i_nim=2,size(images)
       tmp = 0d0
       do i=1,size(images(1)%M%modes_v,2)
          l1 = calc_ang(images(i_nim-1)%M%modes_v(:,i),images(i_nim)%M%modes_v(:,i))
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
! I have removed the first configurations spini and spinf but we can put it back (???)

subroutine geodesic_path(amp_rnd,images)
    use m_get_random
    implicit none
    real(8), intent(in)             :: amp_rnd
    type(lattice), intent(inout)    :: images(:)
!   real(8), target, intent(inout) :: path(:,:,:)
    ! internal variable
    integer             :: N_cell,nim,dim_mode_mag
    real(8)             :: norm_local,dumy(3)
    integer             :: i,i_nim,j
    real(8), allocatable, dimension(:,:)    :: path_one
    real(8), allocatable, dimension(:)      :: u(:),ax(:)
    
    N_cell=images(1)%Ncell
    nim=size(images)

    dim_mode_mag=images(1)%M%dim_mode
    
    allocate(path_one(dim_mode_mag,nim),u(dim_mode_mag),ax(dim_mode_mag),source=0.0d0)
    
    !could certainly be optimized with more vector instructions moving parts up from low-laying loops 
    do i_nim = 1,nim
        do i=1,N_cell
            call geodesic_path_one(nim,images(1)%M%modes_v(:,i),images(nim)%M%modes_v(:,i),ax,path_one,dim_mode_mag)
            do j=1,dim_mode_mag/3
                u(1)=2d0*(get_rand_classic()-1d0)
                u(2)=2d0*(get_rand_classic()-1d0)
                u(3)=2d0*(get_rand_classic()-1d0)
                dumy=path_one(3*(j-1)+1:3*j ,i_nim)+amp_rnd*u
                norm_local=norm(dumy)
                if (norm_local.gt.1.0d-8) images(i_nim)%M%modes_v(3*(j-1)+1:3*j,i)=dumy/norm_local
            enddo
        end do
    end do

end subroutine geodesic_path

end module m_path
