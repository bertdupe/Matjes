module m_propagator
use m_constants, only : hbar
use m_vector, only : cross,norm
implicit none

public

contains

! ----------------------------------------------
! ----------------------------------------------
! LLG equation of motion that returns the D_M
function LLG(B,damping,spini,size_B)
use m_torques, only : update_B
implicit none
integer, intent(in) :: size_B
real(kind=8), intent(in) :: spini(:),damping,B(:)   !,stm_field_torque
real(kind=8) :: LLG(size_B)
!dummy
real(kind=8) :: stepdamp(size_B),norm_S,S_norm(size_B)
real(kind=8) :: LLG_int(size_B)

LLG=0.0d0
LLG_int=0.0d0

norm_S=norm(spini)
S_norm=spini/norm_S

stepdamp=cross(S_norm,B,1,size_b)

LLG_int=-B-damping*stepdamp

call update_B(S_norm,damping,LLG_int)

LLG=cross(S_norm,LLG_int,1,size_b)

end function LLG

! ----------------------------------------------
! ----------------------------------------------
! LLG equation of motion that returns the B used in MxB
function LLG_B(B,damping,spini,size_B)
use m_torques, only : update_B
implicit none
integer, intent(in) :: size_B
real(kind=8), intent(in) :: spini(:),damping,B(:)   !,stm_field_torque
real(kind=8) :: LLG_B(size_B)
!dummy
real(kind=8) :: stepdamp(size_B),norm_S,S_norm(size_B)

LLG_B=0.0d0

norm_S=norm(spini)
S_norm=spini/norm_S

stepdamp=cross(S_norm,B,1,size_b)

LLG_B=-B-damping*stepdamp

call update_B(S_norm,damping,LLG_B)

end function LLG_B

end module m_propagator
