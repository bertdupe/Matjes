       module m_solver
        interface integrate
         module procedure simple
         module procedure integrate_pred
         module procedure integrate_SIB
         module procedure integrate_SIB_NC_ohneT
        end interface integrate
        interface minimization
         module procedure euler_minimization
         module procedure euler_o2_minimization
        end interface minimization
       contains

! ----------------------------------------------
! SIA prediction integration scheme. Not norm conserving
       function integrate_pred(timestep,spin1,B,kt,damping,stmtemp,state, &
      & i_torque,stmtorque,torque_FL,torque_AFL,adia,nonadia,storque,maxh,Ipol,i_x,i_y,i_z,i_m,spin)
       use m_constants, only : hbar
       use m_randist
       use m_dynamic, only : htor
       use m_lattice, only : tableNN,masque
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
       type(mtprng_state), intent(inout) :: state
       ! dummy
       real(kind=8) :: Dtemp,W(3),dt
       real(kind=8) :: step(3),steptor(3),stepadia(3),stepsttor(3),steptemp(3),sum_step(3)
       integer :: v_x,v_y,v_z,v_m
#ifndef CPP_MPI
       integer, parameter :: Xstart=1
       integer, parameter :: Ystart=1
       integer, parameter :: Zstart=1
#endif

       if (masque(1,i_x,i_y,i_z).eq.0) then
         integrate_pred=0.0d0
         return
       endif

       step=0.0d0
       steptor=0.0d0
       stepadia=0.0d0
       stepsttor=0.0d0
       steptemp=0.0d0

       dt=timestep/hbar/(1+damping**2)

        step=cross(B,spin1)
        if (i_torque) steptor=cross(spin1,Ipol)
        if (i_torque) then
         v_x=tableNN(1,1,i_x,i_y,i_z,i_m)
         v_y=tableNN(2,1,i_x,i_y,i_z,i_m)
         v_z=tableNN(3,1,i_x,i_y,i_z,i_m)
         v_m=tableNN(4,1,i_x,i_y,i_z,i_m)
         stepadia=cross(-spin(4:6,v_x,v_y,v_z,v_m),spin1)
        endif
        if (stmtorque) stepsttor=cross(spin1,Ipol*htor(i_x,i_y,i_z))

       if (kt.gt.1.0d-10) then
        Dtemp=damping/(1+damping**2)*kT/sqrt(B(1)**2+B(2)**2+B(3)**2)
        if (stmtemp) then
         W=(/randist(1,state),randist(1,state),randist(1,state)/)
         steptemp=(sqrt(2.0d0*Dtemp)*W+damping*sqrt(2.0d0*Dtemp)*cross(spin(4:6,i_x,i_y,i_z,i_m),W))*htor(i_x,i_y,i_z)/maxh
        else
         W=(/randist(1,state),randist(1,state),randist(1,state)/)
         steptemp=+sqrt(2.0d0*Dtemp)*W+damping*sqrt(2.0d0*Dtemp)*cross(spin(4:6,i_x,i_y,i_z,i_m),W)
        endif
       endif

       sum_step=step+damping*cross(spin(4:6,i_x,i_y,i_z,i_m),step)+                        &
     &     torque_FL*(1.0d0-damping*torque_AFL)*steptor+(damping+torque_AFL)*torque_FL*    &
     &     cross(spin(4:6,i_x,i_y,i_z,i_m),steptor)+adia*                                  &
     &     cross(spin(4:6,i_x,i_y,i_z,i_m),stepadia)-nonadia*stepadia                      &
     &     +storque*cross(stepsttor,spin(4:6,i_x,i_y,i_z,i_m))

       step=spin1+(sum_step*dt+cross(spin1,steptemp)*sqrt(dt))/2.0d0

       integrate_pred=step/sqrt(step(1)**2+step(2)**2+step(3)**2)

       end function integrate_pred

! ----------------------------------------------
! SIB integration scheme. the 3x3 is inverted by hand
       function integrate_SIB(timestep,spin1,B,kt,damping,stmtemp,state, &
      & i_torque,stmtorque,torque_FL,torque_AFL,adia,nonadia,storque,maxh,check,Ipol,i_x,i_y,i_z,i_m,spin)
      use m_constants, only : hbar
       use m_randist
       use m_dynamic, only : htor
       use m_lattice, only : masque,tableNN
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
       type(mtprng_state), intent(inout) :: state
       ! dummy
       real(kind=8) :: Dtemp,W(3),droite(3),dt,denominator,dumy
       real(kind=8) :: step(3),steptor(3),stepadia(3),stepsttor(3),sum_step(3),steptemp(3)
       integer :: v_x,v_y,v_z,v_m
#ifndef CPP_MPI
       integer, parameter :: Xstart=1
       integer, parameter :: Ystart=1
       integer, parameter :: Zstart=1
#endif

       if (masque(1,i_x,i_y,i_z).eq.0) then
         integrate_SIB=0.0d0
         return
       endif

       step=0.0d0
       steptor=0.0d0
       stepadia=0.0d0
       stepsttor=0.0d0
       w=0.0d0
       steptemp=0.0d0

        step=cross(B,spin1)
        if (i_torque) steptor=cross(spin1,Ipol)
        if (i_torque) then
         v_x=tableNN(1,1,i_x,i_y,i_z,i_m)
         v_y=tableNN(2,1,i_x,i_y,i_z,i_m)
         v_z=tableNN(3,1,i_x,i_y,i_z,i_m)
         v_m=tableNN(4,1,i_x,i_y,i_z,i_m)
         stepadia=cross(-spin(4:6,v_x,v_y,v_z,v_m),spin1)
        endif
        if (stmtorque) stepsttor=cross(spin1,Ipol*htor(i_x,i_y,i_z))

       if (kt.gt.1.0d-10) then
        Dtemp=damping/(1+damping**2)*kT/sqrt(B(1)**2+B(2)**2+B(3)**2)
        if (stmtemp) then
         W=(/randist(1,state),randist(1,state),randist(1,state)/)
         steptemp=(sqrt(2.0d0*Dtemp)*W+damping*sqrt(2.0d0*Dtemp)*cross(spin(4:6,i_x,i_y,i_z,i_m),W))*htor(i_x,i_y,i_z)/maxh
        else
         W=(/randist(1,state),randist(1,state),randist(1,state)/)
         steptemp=sqrt(2.0d0*Dtemp)*W+damping*sqrt(2.0d0*Dtemp)*cross(spin(4:6,i_x,i_y,i_z,i_m),W)
        endif
       endif

       sum_step=step+damping*cross(spin(4:6,i_x,i_y,i_z,i_m),step)+     &
     &     torque_FL*(1.0d0-damping*torque_AFL)*steptor+                &
     &     torque_FL*(damping+torque_AFL)*                              &
     &     cross(spin(4:6,i_x,i_y,i_z,i_m),steptor)+adia*               &
     &     cross(spin(4:6,i_x,i_y,i_z,i_m),stepadia)-nonadia*stepadia   &
     &     +storque*cross(stepsttor,spin(4:6,i_x,i_y,i_z,i_m))

       dt=timestep/hbar/(1+damping**2)


! 3x3 system
! SX=droite
       droite=spin1+sum_step*dt/2.0d0+sqrt(dt)/2.0d0*cross(W,steptemp)
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
       check(1)=check(1)+sum(cross(integrate_SIB,B)**2)
       check(2)=check(2)+dot_product(integrate_SIB,B)
!!! end check

       end function integrate_SIB

! ----------------------------------------------
! SIB integration scheme for predicator. norm conserving
       function integrate_SIB_NC_ohneT(timestep,spin1,B,damping, &
      & i_torque,stmtorque,torque_FL,torque_AFL,adia,nonadia,storque,maxh,Ipol,i_x,i_y,i_z,i_m,spin)
       use m_constants, only : hbar
       use m_randist
       use m_dynamic, only : htor
       use m_lattice, only : tableNN,masque
       use m_vector, only : cross,norm
       use mtprng
#ifdef CPP_MPI
       use m_make_box, only : Xstart,Xstop,Ystart,Ystop,Zstart,Zstop
#endif
       implicit none
       real(kind=8) :: integrate_SIB_NC_ohneT(3)
       real(kind=8), intent(in) :: timestep,B(:),spin1(:),damping,torque_FL,torque_AFL,adia, &
      & nonadia,storque,maxh,Ipol(:),spin(:,:,:,:,:)
       integer, intent(in) :: i_x,i_y,i_z,i_m
       logical, intent(in) :: i_torque,stmtorque
       ! dummy
       real(kind=8) :: sum_step(3),droite(3),denominator,dumy
       real(kind=8) :: step(3),steptor(3),stepadia(3),stepsttor(3),dt
       integer :: v_x,v_y,v_z,v_m
#ifndef CPP_MPI
       integer, parameter :: Xstart=1
       integer, parameter :: Ystart=1
       integer, parameter :: Zstart=1
#endif

       if (masque(1,i_x,i_y,i_z).eq.0) then
         integrate_SIB_NC_ohneT=0.0d0
         return
       endif

       dt=timestep/hbar/(1+damping**2)

       step=0.0d0
       steptor=0.0d0
       stepadia=0.0d0
       stepsttor=0.0d0

        step=cross(B,spin1)
        if (i_torque) steptor=cross(spin1,Ipol)
        if (i_torque) then
         v_x=tableNN(1,1,i_x,i_y,i_z,i_m)
         v_y=tableNN(2,1,i_x,i_y,i_z,i_m)
         v_z=tableNN(3,1,i_x,i_y,i_z,i_m)
         v_m=tableNN(4,1,i_x,i_y,i_z,i_m)
         stepadia=cross(-spin(4:6,v_x,v_y,v_z,v_m),spin1)
        endif
        if (stmtorque) stepsttor=cross(spin1,Ipol*htor(i_x,i_y,i_z))

       sum_step=step+damping*cross(spin(4:6,i_x,i_y,i_z,i_m),step)+     &
     &     torque_FL*(1.0d0-damping*torque_AFL)*steptor+                                           &
     &     torque_FL*(damping+torque_AFL)*                                        &
     &     cross(spin(4:6,i_x,i_y,i_z,i_m),steptor)+adia*               &
     &     cross(spin(4:6,i_x,i_y,i_z,i_m),stepadia)-nonadia*stepadia   &
     &     +storque*cross(stepsttor,spin(4:6,i_x,i_y,i_z,i_m))

! 3x3 system
! SX=droite
       droite=spin1+sum_step*dt/2.0d0
       denominator=(4.0d0+B(1)**2*dt**2+B(2)**2*dt**2+B(3)**2*dt**2)

! first term
       integrate_SIB_NC_ohneT(1)=-(-4.0d0*droite(1)+ &
     &  dt*(2.0d0*B(3)*droite(2)-2.0d0*B(2)*droite(3)) + &
     &  dt**2*(-B(1)**2*droite(1)-B(1)*B(2)*droite(2)-B(1)*B(3)*droite(3))) &
     &  /denominator

! second term
       integrate_SIB_NC_ohneT(2)=-(-4.0d0*droite(2)+ &
     &  dt*(-2.0d0*B(3)*droite(1)+2.0d0*B(1)*droite(3))+ &
     &  dt**2*(-B(2)**2*droite(2)-B(1)*B(2)*droite(1)-B(2)*B(3)*droite(3))) &
     &  /denominator

! third term
       integrate_SIB_NC_ohneT(3)=-(-4.0d0*droite(3)+ &
     &  dt*(2.0d0*B(2)*droite(1)-2.0d0*B(1)*droite(2))+ &
     &  dt**2*(-B(3)**2*droite(3)-B(1)*B(3)*droite(1)-B(2)*B(3)*droite(2))) &
     &  /denominator


       end function integrate_SIB_NC_ohneT


! ----------------------------------------------
! ----------------------------------------------
! Euler integration scheme
       function simple(timestep,B,kt,stmtemp,maxh,i_x,i_y,i_z,i_m,damping,Ipol,torque_FL,torque_AFL,adia,nonadia, &
      & storque,i_torque,stmtorque,spin,shape_spin,tableNN)
       use m_constants, only : hbar
       use m_randist
       use m_vector, only : cross
       use m_dynamic, only : htor
       implicit none
       real(kind=8) :: simple(3)
       integer, intent(in) :: tableNN(:,:,:,:,:,:),shape_spin(:)
       real(kind=8), intent(in) :: timestep,B(:),kt,maxh,spin(:,:,:,:,:)
       real(kind=8), intent(in) :: damping,torque_FL,torque_AFL,adia,nonadia,storque
       real(kind=8), intent(in) :: Ipol(:)
       integer, intent(in) :: i_x,i_y,i_z,i_m
       logical, intent(in) :: stmtemp,i_torque,stmtorque
       !dummy
       real(kind=8) :: spinfin(3),steptemp(3),step(3),steptor(3),stepadia(3),stepsttor(3),stepdamp(3),spini(3)
       real(kind=8) :: dt,ds(3),norm_S,S_norm(3)
       integer :: v_x,v_y,v_z,v_m
       integer :: X,Y,Z

       steptemp=0.0d0
       dt=timestep/hbar/(1+damping**2)
       simple=0.0d0
       step=0.0d0
       steptor=0.0d0
       stepadia=0.0d0
       stepsttor=0.0d0
       X=shape_spin(1)-2
       Y=shape_spin(1)-1
       Z=shape_spin(1)-0
       spini=spin(X:Z,i_x,i_y,i_z,i_m)
       stepdamp=0.0d0

       if (kt.gt.1.0d-10) then
        if (stmtemp) then
         steptemp=cross(spini,(/randist(kt),randist(kt),randist(kt)/))*htor(i_x,i_y,i_z)/maxh
        else
         steptemp=cross(spini,(/randist(kt),randist(kt),randist(kt)/))
        endif
       endif

       v_x=tableNN(1,1,i_x,i_y,i_z,i_m)
       v_y=tableNN(2,1,i_x,i_y,i_z,i_m)
       v_z=tableNN(3,1,i_x,i_y,i_z,i_m)
       v_m=tableNN(4,1,i_x,i_y,i_z,i_m)

       norm_S=sqrt(spini(1)**2+spini(2)**2+spini(3)**2)
       S_norm=spini/norm_S


       if (i_torque) steptor=cross(S_norm,Ipol)
       if (i_torque) stepadia=cross(-spin(:,v_x,v_y,v_z,v_m),S_norm)
       if (stmtorque) stepsttor=cross(S_norm,Ipol*htor(i_x,i_y,i_z))

       step=cross(B,S_norm)

       stepdamp=cross(S_norm,step)

       ds=step+damping*stepdamp

       ds=ds+torque_FL*(1.0d0-damping*torque_AFL)*steptor+     &
     &   torque_FL*(torque_AFL+damping)*cross(S_norm,steptor)+adia*       &
     &   cross(S_norm,stepadia)-nonadia*stepadia+storque*cross(stepsttor,S_norm)+steptemp


       spinfin=S_norm+dt*ds

       simple=spinfin/sqrt(spinfin(1)**2+spinfin(2)**2+spinfin(3)**2)

       end function simple

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
