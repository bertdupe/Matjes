module m_solver
use m_derived_types

   interface minimization
      module procedure euler_minimization
      module procedure euler_o2_minimization
   end interface minimization

private
public :: Euler,euler_old,implicite,minimization
contains

!
! dt comes in in units of fs
!
!
!

! ----------------------------------------------
! ----------------------------------------------
! implicit integration scheme
function implicite(mode_t,D_mode,DT_mode,dt,size_mode)
use m_constants, only : hbar
use m_vector, only : cross,norm
implicit none
integer, intent(in) :: size_mode
real(kind=8), intent(in) :: mode_t(:),D_mode(:),dt,DT_mode(:)
real(kind=8) :: implicite(size_mode)
!dummy
real(kind=8) :: droite(size_mode),denominator,B_int(size_mode)
integer :: i,start

implicite=0.0d0
B_int=D_mode*dt/hbar

do i=1,size_mode/3
  start=3*(i-1)+1

! 3x3 system
! SX=droite

  droite(start:start+2)=mode_t(start:start+2)+cross(mode_t,D_mode,start,start+2)*dt/hbar
  denominator=(1.0d0+norm(B_int(start:start+2))**2)

! first term
  implicite(start)=-(-droite(start)+ &
     &  (B_int(start+2)*droite(start+1)-B_int(start+1)*droite(start+2)) + &
     &  (-B_int(start)**2*droite(start)-B_int(start)*B_int(start+1)*droite(start+1)-B_int(start)*B_int(start+2)*droite(start+2))) &
     &  /denominator

! second term
  implicite(start+1)=-(-droite(start+1)+ &
     &  (-B_int(start+2)*droite(start)+B_int(start)*droite(start+2))+ &
     &  (-B_int(start+1)**2*droite(start+1)-B_int(start)*B_int(start+1)*droite(start)-B_int(start+1)*B_int(start+2)*droite(start+2))) &
     &  /denominator

! third term
  implicite(start+2)=-(-droite(start+2)+ &
     &  (B_int(start+1)*droite(start)-B_int(start)*droite(start+1))+ &
     &  (-B_int(start+2)**2*droite(start+2)-B_int(start)*B_int(start+2)*droite(start)-B_int(start+1)*B_int(start+2)*droite(start+1))) &
     &  /denominator

enddo

end function implicite

! ----------------------------------------------
! ----------------------------------------------
! Euler integration scheme
function euler_old(mode_t,D_mode,DT_mode,dt,size_mode)result(euler)
use m_vector, only : cross,norm
use m_constants, only : hbar
implicit none
integer, intent(in) :: size_mode
real(kind=8), intent(in) :: mode_t(:),D_mode(:),dt,DT_mode(:)
real(kind=8) :: euler(size_mode)
!dummy
real(kind=8) :: euler_int(size_mode)
real(kind=8) :: norm_mode,norm_int
integer :: i,start,end

euler=mode_t

do i=1,size_mode/3
   start=3*(i-1)+1
   end=3*i

   norm_mode=norm(mode_t(start:end))
   euler_int(start:end)=mode_t(start:end)+(D_mode(start:end)*dt+sqrt(dt)*DT_mode(start:end))/hbar
   norm_int=norm(euler_int(start:end))
   if (norm(euler_int(start:end)).gt.1.0d-8) euler(start:end)=euler_int(start:end)*norm_mode/norm_int
enddo

end function 


function euler(m,Dmag_int,DTmag_int,dt)result(Mout)
    !do the LLG step from data provided in large contiguous arrays (with shape (3,Nsite*Nmag))
    use m_constants, only : hbar
    real(8),intent(in),contiguous   ::  M(:,:),Dmag_int(:,:),DTmag_int(:,:)
    real(8),intent(in)              ::  dt
    real(8),allocatable             ::  Mout(:,:)
    !internal
    real(8) :: m_norm(size(M,2))
    real(8) :: int_norm(size(M,2))
    real(8) :: euler_tmp(size(M,1),size(M,2))
    integer :: i

    allocate(Mout(size(M,1),size(M,2)),source=0.0d0)
    !ADD DT_mode if necessary
    if(size(M,1)/=3.or.size(Dmag_int,1)/=3.or.size(DTmag_int,1)/=3.or.size(Mout,1)/=3) ERROR STOP "euler INPUT NEEDS TO BE 3-vector"
    m_norm=norm2(m,dim=1)
    euler_tmp=M+Dmag_int*dt/hbar+DTmag_int*dt/hbar
    int_norm=norm2(euler_tmp,dim=1)
    mout=M
    forall(i=1:size(int_norm), int_norm(i)>1.0d-8 ) mout(:,i)=euler_tmp(:,i)*m_norm(i)/int_norm(i)

end function




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!
!!!!!!!!!!!!!   Minimization part
!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine euler_minimization(initial,force,predicator,dt,masse)
implicit none
real(kind=8), intent(in) :: force(:),initial(:)
real(kind=8), intent(out) :: predicator(size(initial))
real(kind=8), intent(in) :: dt,masse
! dummy variable

predicator=force*dt/masse+initial

end subroutine euler_minimization
!velocity(:,iomp),(force(:,iomp)+F_temp)/2.0d0,V_ef
subroutine euler_o2_minimization(spin,v_spin,force,predicator,dt,masse)
implicit none
real(kind=8), intent(in) :: force(:),v_spin(:),spin(:)
real(kind=8), intent(in) :: dt,masse
real(kind=8), intent(out) :: predicator(size(spin))
! dummy variable
real(kind=8) :: s_dumy(size(spin)),norm_fin,norm_ini
integer :: i

s_dumy=0.0d0
predicator=0.0d0

s_dumy=force*dt**2/masse/2.0d0+v_spin*dt+spin

do i=1,size(spin)/3
  norm_ini=sqrt(sum(spin(3*(i-1)+1:3*i)**2))
  norm_fin=sqrt(sum(s_dumy(3*(i-1)+1:3*i)**2))
  if (norm_fin.gt.1.0d-8) predicator(3*(i-1)+1:3*i)=s_dumy(3*(i-1)+1:3*i)/norm_fin*norm_ini
enddo

end subroutine euler_o2_minimization

end module m_solver
