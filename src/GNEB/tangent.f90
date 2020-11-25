module m_tangent
use m_basic_types, only : vec
use m_vector, only : norm, project
use m_projection
use m_derived_types, only : lattice
implicit none

interface tang
   module procedure tang_oneimage,tang_spec
   module procedure lattice_single 
   module procedure lattice_special
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

N_cell=size(coo,2)
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

!project tangent onto tangent space and normalize it
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

tau=tau/sqrt(tmp)

end subroutine tang_spec

!> Estimate tangent to the path according to the special definition: doi:10.1016/j.cpc.2015.07.001
subroutine lattice_special(im,images,u,tau)
    integer, intent(in)             :: im
    type(lattice), intent(inout)    :: images(:)
    real(8), intent(in)             :: u(:)  !energies
    real(8), intent(out)            :: tau(:,:)
    ! internal variables
    real(8) :: u1, u2, u3,dumin,dumax,tmp
    real(8), allocatable :: taup(:,:),taum(:,:)
    
    u1=u(im-1)
    u2=u(im)
    u3=u(im+1)
    if (u3>u2.and.u2>u1) then
        tau=images(im+1)%M%modes_v-images(im)%M%modes_v
    elseif (u1>u2.and.u2>u3) then
        tau=images(im)%M%modes_v-images(im-1)%M%modes_v
    else
        allocate(taup,taum,mold=tau)
        taup=images(im+1)%M%modes_v-images(im)%M%modes_v
        taum=images(im)%M%modes_v-images(im-1)%M%modes_v
        dumax=dabs(u3-u2)
        dumin=dabs(u1-u2)
        if (dumax<dumin) then
           tmp=dumax
           dumax=dumin
           dumin=tmp
        end if
        if (u3>u1) then
            tau=dumax*taup+dumin*taum
        else
            tau=dumin*taup+dumax*taum
        end if
    end if
    Call project(tau,images(im)%M%modes_v)
    tau=tau/norm2(tau)
end subroutine 



!
! Calculate the tangent to the path
!
subroutine tang_oneimage(nim,im,coo,tau)
implicit none
integer, intent(in)         :: im,nim
real(kind=8), intent(in)    :: coo(:,:,:)
real(kind=8), intent(out)   :: tau(:,:)
! internal variables
real(kind=8) :: tmp
integer :: i,N_cell,u1,u2

N_cell=size(coo,2)

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

tau = tau/sqrt(tmp)

end subroutine tang_oneimage


!
! Calculate the tangent to the path
!
subroutine lattice_single(im,images,tau)
    integer, intent(in)             :: im
    type(lattice), intent(inout)    :: images(:)
    real(8), intent(out)            :: tau(:,:)
    ! internal variables
    integer :: u1,u2
    
    u1=im+1
    u2=im-1
    if (im==1) u2=im
    if (im==size(images)) u1=im
    
    tau=images(u1)%M%modes_v-images(u2)%M%modes_v
    Call project(tau,images(im)%M%modes_v)
    tau=tau/norm2(tau)
end subroutine 


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
