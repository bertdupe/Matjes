module m_solver
use m_derived_types

   interface minimization
      module procedure euler_minimization
      module procedure euler_o2_minimization
   end interface minimization

private
public :: Euler,implicite,minimization
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

write(*,*) mode_t(:)
write(*,*) D_mode
write(*,*) DT_mode

do i=1,size_mode/3
  start=3*(i-1)+1

! 3x3 system
! SX=droite

  droite(start:start+2)=mode_t(start:start+2)+cross(mode_t,D_mode,start,start+2)*dt/hbar
  denominator=(1.0d0+norm(B_int(start:start+2))**2)

! first term
  implicite(start)=-(-droite(start)+ &
     &  (-B_int(start+2)*droite(start+1)+B_int(start+1)*droite(start+2)) + &
     &  (-B_int(start)**2*droite(start)-B_int(start)*B_int(start+1)*droite(start+1)-B_int(start)*B_int(start+2)*droite(start+2))) &
     &  /denominator

! second term
  implicite(start+1)=-(-droite(start+1)+ &
     &  (B_int(start+2)*droite(start)-B_int(start)*droite(start+2))+ &
     &  (-B_int(start+1)**2*droite(start+1)-B_int(start)*B_int(start+1)*droite(start)-B_int(start+1)*B_int(start+2)*droite(start+2))) &
     &  /denominator

! third term
  implicite(start+2)=-(-droite(start+2)+ &
     &  (-B_int(start+1)*droite(start)+B_int(start)*droite(start+1))+ &
     &  (-B_int(start+2)**2*droite(start+2)-B_int(start)*B_int(start+2)*droite(start)-B_int(start+1)*B_int(start+2)*droite(start+1))) &
     &  /denominator

enddo
write(*,*) 'spin', implicite(start:start+2),sqrt(sum(implicite(start:start+2)**2))


end function implicite

! ----------------------------------------------
! ----------------------------------------------
! Euler integration scheme
function euler(mode_t,D_mode,DT_mode,dt,size_mode)
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

!if (BT_norm.gt.1.0d-10) then
!   step_T=cross(S_norm,BT,1,3)
!   DTs= step_T+damping*cross(S_norm,step_T,1,3)
!endif

do i=1,size_mode/3
   start=3*(i-1)+1
   end=3*i
   norm_mode=norm(mode_t(start:end))
   euler_int(start:end)=mode_t(start:end)+(D_mode(start:end)*dt+sqrt(dt)*DT_mode(start:end))/hbar
   norm_int=norm(euler_int(start:end))
   if (norm(euler_int(start:end)).gt.1.0d-8) euler(start:end)=euler_int(start:end)*norm_mode/norm_int
enddo

end function euler


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!
!!!!!!!!!!!!!   Minimization part
!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine euler_minimization(initial,force,predicator,dt,masse,N)
implicit none
integer, intent(in) :: N
real(kind=8), intent(out) :: predicator(N)
real(kind=8), intent(in) :: force(N),initial(N)
real(kind=8), intent(in) :: dt,masse
! dummy variable

predicator=force*dt/masse+initial

end subroutine euler_minimization

subroutine euler_o2_minimization(spin,v_spin,force,predicator,dt,masse,N)
implicit none
integer, intent(in) :: N
real(kind=8), intent(out) :: predicator(N)
real(kind=8), intent(in) :: force(N),v_spin(N),spin(N)
real(kind=8), intent(in) :: dt,masse
! dummy variable
real(kind=8) :: s_dumy(N),norm

s_dumy=0.0d0
predicator=0.0d0

s_dumy=force*dt**2/masse/2.0d0+v_spin*dt+spin
norm=sqrt(sum(s_dumy**2))

predicator=s_dumy/norm

end subroutine euler_o2_minimization

end module m_solver
