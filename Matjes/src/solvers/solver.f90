module m_solver
use m_derived_types

   interface minimization
      module procedure euler_minimization
      module procedure euler_o2_minimization
   end interface minimization

private
public :: Euler,integrate_SIB_NC_ohneT,implicite,minimization
contains

!
! dt comes in in units of fs
!
!
!

! ----------------------------------------------
! SIA prediction integration scheme. Not norm conserving
       function integrate_pred(timestep,spin1,B,kt,damping,stmtemp, &
      & i_torque,stmtorque,torque_FL,torque_AFL,adia,nonadia,storque,maxh,Ipol,i_x,i_y,i_z,i_m,spin)
       use m_constants, only : hbar
       use m_randist
       use m_vector, only : cross,norm
       use mtprng
#ifdef CPP_MPI
       use m_make_box, only : Xstart,Xstop,Ystart,Ystop,Zstart,Zstop
#endif
       implicit none
       real(kind=8) :: integrate_pred(3)
       real(kind=8), intent(in) :: timestep,B(:),spin1(:),damping,torque_FL,torque_AFL,adia, &
      & nonadia,storque,maxh,Ipol(:),kt,spin(:,:,:,:,:)
       integer, intent(in) :: i_x,i_y,i_z,i_m
       logical, intent(in) :: stmtemp,i_torque,stmtorque
       ! dummy
       type(mtprng_state) :: state
       real(kind=8) :: Dtemp,W(3),dt
       real(kind=8) :: step(3),steptor(3),stepadia(3),stepsttor(3),steptemp(3),sum_step(3)
       integer :: v_x,v_y,v_z,v_m
#ifndef CPP_MPI
       integer, parameter :: Xstart=1
       integer, parameter :: Ystart=1
       integer, parameter :: Zstart=1
#endif

       step=0.0d0
       steptor=0.0d0
       stepadia=0.0d0
       stepsttor=0.0d0
       steptemp=0.0d0

       dt=timestep/hbar/(1+damping**2)

        step=cross(B,spin1,1,3)
        if (i_torque) steptor=cross(spin1,Ipol,1,3)

       if (kt.gt.1.0d-10) then
        Dtemp=damping/(1+damping**2)*kT/norm(B)
        if (stmtemp) then
!         W=(/randist(state),randist(state),randist(state)/)
!         steptemp=(sqrt(2.0d0*Dtemp)*W+damping*sqrt(2.0d0*Dtemp)*cross(spin(:,i_x,i_y,i_z,i_m),W,1,3))*htor(i_x,i_y,i_z)/maxh
        else
         W=(/randist(state),randist(state),randist(state)/)
         steptemp=+sqrt(2.0d0*Dtemp)*W+damping*sqrt(2.0d0*Dtemp)*cross(spin(:,i_x,i_y,i_z,i_m),W,1,3)
        endif
       endif

       sum_step=step+damping*cross(spin(:,i_x,i_y,i_z,i_m),step,1,3)+                        &
     &     torque_FL*(1.0d0-damping*torque_AFL)*steptor+(damping+torque_AFL)*torque_FL*    &
     &     cross(spin(:,i_x,i_y,i_z,i_m),steptor,1,3)+adia*                                  &
     &     cross(spin(:,i_x,i_y,i_z,i_m),stepadia,1,3)-nonadia*stepadia                      &
     &     +storque*cross(stepsttor,spin(:,i_x,i_y,i_z,i_m),1,3)

       step=spin1+(sum_step*dt+cross(spin1,steptemp,1,3)*sqrt(dt))/2.0d0

       integrate_pred=step/norm(step)

       end function integrate_pred

! ----------------------------------------------
! SIB integration scheme. the 3x3 is inverted by hand
       function integrate_SIB(timestep,spin1,B,kt,damping,stmtemp, &
      & i_torque,stmtorque,torque_FL,torque_AFL,adia,nonadia,storque,maxh,check,Ipol,i_x,i_y,i_z,i_m,spin)
      use m_constants, only : hbar
       use m_randist
       use m_vector, only : cross,norm
       use mtprng
#ifdef CPP_MPI
       use m_make_box, only : Xstart,Xstop,Ystart,Ystop,Zstart,Zstop
#endif
       implicit none
       real(kind=8) :: integrate_SIB(3)
       real(kind=8), intent(in) :: timestep,B(:),spin1(:),damping,maxh,kt,torque_FL,torque_AFL,adia,nonadia,storque,Ipol(:),spin(:,:,:,:,:)
       integer, intent(in) :: i_x,i_y,i_z,i_m
       logical, intent(in) :: stmtemp,i_torque,stmtorque
       real(kind=8), intent(inout) :: check(:)
       ! dummy
       type(mtprng_state) :: state
       real(kind=8) :: Dtemp,W(3),droite(3),dt,denominator
       real(kind=8) :: step(3),steptor(3),stepadia(3),stepsttor(3),sum_step(3),steptemp(3)
       integer :: v_x,v_y,v_z,v_m
#ifndef CPP_MPI
       integer, parameter :: Xstart=1
       integer, parameter :: Ystart=1
       integer, parameter :: Zstart=1
#endif

       step=0.0d0
       steptor=0.0d0
       stepadia=0.0d0
       stepsttor=0.0d0
       w=0.0d0
       steptemp=0.0d0

        step=cross(B,spin1,1,3)
        if (i_torque) steptor=cross(spin1,Ipol,1,3)

       if (kt.gt.1.0d-10) then
        Dtemp=damping/(1+damping**2)*kT/norm(B)
        if (stmtemp) then
!         W=(/randist(state),randist(state),randist(state)/)
!         steptemp=(sqrt(2.0d0*Dtemp)*W+damping*sqrt(2.0d0*Dtemp)*cross(spin(:,i_x,i_y,i_z,i_m),W,1,3))*htor(i_x,i_y,i_z)/maxh
        else
         W=(/randist(state),randist(state),randist(state)/)
         steptemp=sqrt(2.0d0*Dtemp)*W+damping*sqrt(2.0d0*Dtemp)*cross(spin(:,i_x,i_y,i_z,i_m),W,1,3)
        endif
       endif

       sum_step=step+damping*cross(spin(:,i_x,i_y,i_z,i_m),step,1,3)+     &
     &     torque_FL*(1.0d0-damping*torque_AFL)*steptor+                &
     &     torque_FL*(damping+torque_AFL)*                              &
     &     cross(spin(:,i_x,i_y,i_z,i_m),steptor,1,3)+adia*               &
     &     cross(spin(:,i_x,i_y,i_z,i_m),stepadia,1,3)-nonadia*stepadia   &
     &     +storque*cross(stepsttor,spin(:,i_x,i_y,i_z,i_m),1,3)

       dt=timestep/hbar/(1+damping**2)


! 3x3 system
! SX=droite
       droite=spin1+sum_step*dt/2.0d0+sqrt(dt)/2.0d0*cross(W,steptemp,1,3)
       denominator=(4.0d0+B(1)**2*dt**2+B(2)**2*dt**2+B(3)**2*dt**2+2.0d0*B(1)*W(1)*sqrt(dt)**3+2.0d0*B(2)*W(2)*sqrt(dt)**3 &
     &  +2.0d0*B(3)*W(3)*sqrt(dt)**3+W(1)**2*dt+W(2)**2*dt+W(3)**2*dt)

! first term
       integrate_SIB(1)=-(-4.0d0*droite(1)+ &
     &  dt*(2.0d0*B(3)*droite(2)-droite(3)*W(3)*W(1)-2.0d0*B(2)*droite(3)-droite(1)*W(1)**2-droite(2)*W(2)*W(1)) + &
     &  dt**2*(-B(1)**2*droite(1)-B(1)*B(2)*droite(2)-B(1)*B(3)*droite(3))+ &
     &  sqrt(dt)*2.0d0*(-droite(3)*W(2)+droite(2)*W(3))+ &
     &  sqrt(dt)**3*(-2.0d0*B(1)*droite(1)*W(1)-B(2)*droite(2)*W(1)-B(3)*droite(3)*W(1)-B(1)*droite(2)*W(2)-B(2)*droite(3)*W(3))) &
     &  /denominator

! second term
       integrate_SIB(2)=-(-4.0d0*droite(2)+ &
     &  dt*(-2.0d0*B(3)*droite(1)+2.0d0*B(1)*droite(3)-droite(2)*W(2)**2-droite(1)*W(2)*W(1)-droite(3)*W(2)*W(3))+ &
     &  dt**2*(-B(2)**2*droite(2)-B(1)*B(2)*droite(1)-B(2)*B(3)*droite(3))+ &
     &  sqrt(dt)*2.0d0*(droite(3)*W(1)-droite(1)*W(3))+ &
     &  sqrt(dt)**3*(-B(2)*droite(1)*W(1)-B(1)*droite(1)*W(2)-2.0d0*B(2)*droite(2)*W(2)-B(1)*droite(2)*W(2)-B(3)*droite(3)*W(2))) &
     &  /denominator

! third term
       integrate_SIB(3)=-(-4.0d0*droite(3)+ &
     &  dt*(2.0d0*B(2)*droite(1)-2.0d0*B(1)*droite(2)-droite(3)*W(3)**2-droite(1)*W(3)*W(1)-droite(2)*W(2)*W(3))+ &
     &  dt**2*(-B(3)**2*droite(3)-B(1)*B(3)*droite(1)-B(2)*B(3)*droite(2))+ &
     &  sqrt(dt)*2.0d0*(-droite(2)*W(1)+droite(1)*W(2))+ &
     &  sqrt(dt)**3*(-B(3)*droite(1)*W(1)-B(3)*droite(2)*W(2)-2.0d0*B(3)*droite(3)*W(3)-B(1)*droite(1)*W(3)-B(2)*droite(2)*W(3))) &
     &  /denominator

! force norm to 1
!       dumy=sqrt(integrate_pred_heun(1)**2+integrate_pred_heun(2)**2+integrate_pred_heun(3)**2)
!       integrate_pred_heun=integrate_pred_heun/dumy

! test part without thermal noise
!       droite=spin1+step*dt/2.0d0
!       denominator=(4.0d0+B(1)**2*dt**2+B(2)**2*dt**2+B(3)**2*dt**2)
!
!       integrate_spin_NC(1)=-(-4.0d0*droite(1)+ &
!     &  dt*(2.0d0*B(3)*droite(2)+2.0d0*B(2)*droite(3)) + &
!     &  dt**2*(-B(1)**2*droite(1)-B(1)*B(2)*droite(2)-B(1)*B(3)*droite(3)) &
!     &  )/denominator
!
!       integrate_spin_NC(2)=-(-4.0d0*droite(2)+ &
!     &  dt*(-2.0d0*B(3)*droite(1)+2.0d0*B(1)*droite(3))+ &
!     &  dt**2*(-B(2)**2*droite(2)-B(1)*B(2)*droite(1)-B(2)*B(3)*droite(3)) &
!     &  )/denominator
!
!       integrate_spin_NC(3)=-(-4.0d0*droite(3)+ &
!     &  dt*(2.0d0*B(2)*droite(1)-2.0d0*B(1)*droite(2))+ &
!     &  dt**2*(-B(3)**2*droite(3)-B(1)*B(3)*droite(1)-B(2)*B(3)*droite(2)) &
!     &  )/denominator
!
!       dumy=sqrt(integrate_spin_NC(1)**2+integrate_spin_NC(2)**2+integrate_spin_NC(3)**2)
!       integrate_spin_NC=integrate_spin_NC/dumy

! the temperature is checked with 1 temperature step before
!!! check temperature
       check(1)=check(1)+sum(cross(integrate_SIB,B,1,3)**2)
       check(2)=check(2)+dot_product(integrate_SIB,B)
!!! end check

       end function integrate_SIB

! ----------------------------------------------
! SIB integration scheme for predicator. norm conserving
function integrate_SIB_NC_ohneT(timestep,B,BT,damping,spin)
use m_constants, only : hbar
use m_vector, only : cross,norm
use m_torques
#ifdef CPP_MPI
use m_make_box, only : Xstart,Xstop,Ystart,Ystop,Zstart,Zstop
#endif
implicit none
real(kind=8) :: integrate_SIB_NC_ohneT(3)
real(kind=8), intent(in) :: timestep,B(:),spin(:),damping,BT(:)
! dummy
real(kind=8) :: droite(3),denominator,B_int(3)
real(kind=8) :: step(3),stepdamp(3),dt,ds(3)

dt=timestep/hbar
B_int=B
step=cross(B,spin,1,3)

stepdamp=cross(spin,step,1,3)

ds=step+damping*stepdamp

call update_DS(spin,damping,ds)

call update_B(spin,damping,B_int)

! 3x3 system
! SX=droite
droite=spin+ds*dt/2.0d0
denominator=(4.0d0+B_int(1)**2*dt**2+B_int(2)**2*dt**2+B_int(3)**2*dt**2)

! first term
integrate_SIB_NC_ohneT(1)=-(-4.0d0*droite(1)+ &
     &  dt*(2.0d0*B_int(3)*droite(2)-2.0d0*B_int(2)*droite(3)) + &
     &  dt**2*(-B_int(1)**2*droite(1)-B_int(1)*B(2)*droite(2)-B_int(1)*B(3)*droite(3))) &
     &  /denominator

! second term
integrate_SIB_NC_ohneT(2)=-(-4.0d0*droite(2)+ &
     &  dt*(-2.0d0*B_int(3)*droite(1)+2.0d0*B_int(1)*droite(3))+ &
     &  dt**2*(-B_int(2)**2*droite(2)-B_int(1)*B_int(2)*droite(1)-B_int(2)*B_int(3)*droite(3))) &
     &  /denominator

! third term
integrate_SIB_NC_ohneT(3)=-(-4.0d0*droite(3)+ &
     &  dt*(2.0d0*B_int(2)*droite(1)-2.0d0*B_int(1)*droite(2))+ &
     &  dt**2*(-B_int(3)**2*droite(3)-B_int(1)*B_int(3)*droite(1)-B_int(2)*B_int(3)*droite(2))) &
     &  /denominator


end function integrate_SIB_NC_ohneT

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
