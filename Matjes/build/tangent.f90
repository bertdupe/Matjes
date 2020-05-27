module m_tangent
use m_basic_types, only : vec
use m_vector, only : norm
use m_projection

interface tang
   module procedure tang_oneimage,tang_spec
end interface

private
public :: tang

contains

!> Estimate tangent to the path according to the special definition: doi:10.1016/j.cpc.2015.07.001
subroutine tang_spec(im,coo,u,tau)
implicit none
integer, intent(in) :: im
real(kind=8), intent(in) :: coo(:,:,:)
real(kind=8), intent(in) :: u(:)
real(kind=8), intent(out) :: tau(:,:)
! internal variables
real(kind=8) :: u1, u2, u3,dumin,dumax,tmp,tau_tmp(3)
real(kind=8), allocatable :: taup(:,:),taum(:,:)
integer :: i,N_cell

N_cell=size(coo,1)
allocate(taup(3,N_cell),taum(3,N_cell))

u1=u(im-1)
u2=u(im)
u3=u(im+1)
   !print *,'ene:',u1,u2,u3
if (u3>u2.and.u2>u1) then
      !print *,'I am here!!!'
   do i=1,N_cell
      tau(:,i) = coo(:,i,im+1)-coo(:,i,im)
   end do
elseif (u1>u2.and.u2>u3) then
      !print *,'I am here'
   do i=1,N_cell
      tau(:,i) = coo(:,i,im)-coo(:,i,im-1)
   end do
else
   !print *,'I am here!'
   do i=1,N_cell
      taup(:,i) = coo(:,i,im+1)-coo(:,i,im)
      taum(:,i) = coo(:,i,im)-coo(:,i,im-1)
   end do
   dumax=dabs(u3-u2)
   dumin=dabs(u1-u2)
   if (dumax<dumin) then
      tmp=dumax
      dumax=dumin
      dumin=tmp
   end if
   if (u3>u1) then

      call propagate_tau(N_cell,dumax,dumin,taup,taum,tau)

   else
   !print *,'imhere', 'u1 = ', u1, 'u2 = ', u2,'u3 = ', u3
      call propagate_tau(N_cell,dumin,dumax,taup,taum,tau)

   end if
end if

tmp = 0d0
do i=1,N_cell
   tau_tmp(:) = tau(:,i)
   !print *,'tau_tmp:',tau_tmp
   !print *,'norm_vec:',norm_vec(3,tau(:,i_x,i_y,i_z,i_m))
   call project_force(tau_tmp,coo(:,i,im),tau(:,i))
   !print *,'norm_vec:',norm_vec(3,coo(:,i_x,i_y,i_z,i_m,im))
   !print *,' '
   tmp = tmp+norm(tau(:,i))**2

end do

tmp = dsqrt(tmp)
! print *,'tmp:',tmp
do i=1,N_cell
   tau(:,i) = tau(:,i)/tmp
enddo

end subroutine tang_spec

!
! Calculate the tangent to the path
!
subroutine tang_oneimage(nim,im,coo,tau)
implicit none
integer, intent(in) :: im,nim
real(kind=8), intent(in) :: coo(:,:,:)
real(kind=8), intent(out) :: tau(:,:)
! internal variables
real(kind=8) :: tmp
integer :: i,N_cell,u1,u2

N_cell=size(coo,1)

if (im==1) then
   u1=im+1
   u2=im
elseif (im==nim) then
   u1=im
   u2=im-1
else
   u1=im+1
   u2=im-1
end if


tmp = 0d0
do i=1,N_cell
   tau(:,i)=coo(:,i,u1)-coo(:,i,u2)
   call project_force(tau(:,i),coo(:,i,im),tau(:,i))
   tmp = tmp+norm(tau(:,i))**2
end do

tmp = dsqrt(tmp)
do i=1,N_cell
   tau(:,i) = tau(:,i)/tmp
enddo

end subroutine tang_oneimage


!!!!!!!!!!!!!
!
! simple routine to ppropagate tau
!
subroutine propagate_tau(N_cell,dumax,dumin,taup,taum,tau)
implicit none
integer, intent(in) :: N_cell
real(kind=8), intent(in) :: dumax,dumin
real(kind=8), intent(in) :: taup(:,:),taum(:,:)
real(kind=8), intent(out) :: tau(:,:)
! internal variables
integer :: i

do i=1,N_cell
   tau(:,i)=dumax*taup(:,i)+dumin*taum(:,i)
enddo

end subroutine

end module m_tangent
