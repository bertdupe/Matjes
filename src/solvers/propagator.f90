module m_propagator
use m_constants, only : hbar
use m_vector, only : cross,norm
implicit none

private
public :: LLG!, LLG_B,LLG_old

contains

! ----------------------------------------------
! ----------------------------------------------
! LLG equation of motion that returns the D_M
!function LLG_old(B,damping,spini,size_B) result(llg)
!use m_torques, only : update_DS
!implicit none
!integer, intent(in) :: size_B
!real(kind=8), intent(in) :: spini(:),damping,B(:)   !,stm_field_torque
!real(kind=8) :: LLG(size_B)
!!dummy
!real(kind=8) :: stepdamp(size_B),norm_S,S_norm(size_B)
!real(kind=8) :: LLG_int(size_B)
!
!LLG=0.0d0
!LLG_int=0.0d0
!
!norm_S=norm(spini)
!S_norm=spini/norm_S
!
!stepdamp=cross(S_norm,B,1,size_b)
!
!LLG_int=-B-damping*stepdamp
!
!call update_DS(S_norm,damping,LLG_int)
!
!LLG=cross(S_norm,LLG_int,1,size_b)/(1+damping**2)
!
!end function 


subroutine LLG(B,damping,M,Mout)
    !do the LLG step from data provided in large contiguous array (reshaping to (3,Nsite*Nmag))
    real(8),intent(in)                          ::  damping
    real(8),intent(in),target,contiguous        ::  M(:,:),B(:,:)
    real(8),intent(inout),target,contiguous     ::  Mout(:,:)

    !pointers reshaping to (3,:)
    real(8),pointer :: m3(:,:),B3(:,:),m3_out(:,:)
    
    !local data
    real(8),target  ::  M_norm(size(M,1),size(M,2))
    real(8),target  ::  LLG_int(size(M,1),size(M,2))
    real(8),pointer ::  m3_norm(:,:),LLG_int3(:,:)
    integer             :: Nvec

    Nvec=size(M)/3

    !MOVE RESHAPING routine up?
    m3(1:3,1:Nvec)=>M
    m3_norm(1:3,1:Nvec)=>M_norm
    B3(1:3,1:Nvec)=>B
    LLG_int3(1:3,1:Nvec)=>LLG_int
    M3_out(1:3,1:Nvec)=>Mout(:,:)

    Call normalize_M(M3,M3_norm)
    Call cross_new(M3_norm,B3,LLG_int3)
    LLG_int=-B-damping*LLG_int
    Call cross_new(m3_norm,LLG_int3,m3_out)
    Mout(:,:)=Mout(:,:)/(1.0+damping*damping)
    nullify(m3,m3_norm,B3,LLG_int3,M3_out)

    !ADD EXTERNAL TORQUES IF NECESSARY
end subroutine

!help routine for LLG_new
subroutine cross_new(vec1,vec2,vec_out)
    real(8),intent(in)      ::  vec1(:,:)
    real(8),intent(in)      ::  vec2(:,:)
    real(8),intent(out)     ::  vec_out(size(vec1,1),size(vec1,2))

    integer                 ::  i,Nvec
    Nvec=size(vec1,2)
    do i=1,Nvec
        vec_out(1,i)=vec1(2,i)*vec2(3,i)-vec1(3,i)*vec2(2,i)
        vec_out(2,i)=vec1(3,i)*vec2(1,i)-vec1(1,i)*vec2(3,i)
        vec_out(3,i)=vec1(1,i)*vec2(2,i)-vec1(2,i)*vec2(1,i)
    enddo
end subroutine

!help routine for LLG_new
subroutine normalize_M(M,M_norm)
    !normalize vectors
    !first dimension of input has to be vector dimension to be normalized
    
    real(8),intent(in)         ::  M(:,:)
    real(8),intent(inout)      ::  M_norm(:,:)

    real(8)             :: norm(size(M,2))
    integer             :: i

    norm=norm2(m,dim=1)
    do i=1,size(M,2)
        m_norm(:,i)=m(:,i)/norm(i)
    enddo
end subroutine


!! ----------------------------------------------
!! ----------------------------------------------
!! LLG equation of motion that returns the B used in MxB
!function LLG_B(B,damping,spini,size_B)
!use m_torques, only : update_B
!implicit none
!integer, intent(in) :: size_B
!real(kind=8), intent(in) :: spini(:),damping,B(:)   !,stm_field_torque
!real(kind=8) :: LLG_B(size_B)
!!dummy
!real(kind=8) :: stepdamp(size_B),norm_S,S_norm(size_B)
!
!LLG_B=0.0d0
!
!norm_S=norm(spini)
!S_norm=spini/norm_S
!
!stepdamp=cross(S_norm,B,1,size_b)
!
!LLG_B=-B-damping*stepdamp
!
!call update_B(S_norm,damping,LLG_B)
!
!LLG_B=LLG_B/(1+damping**2)
!
!end function LLG_B

end module m_propagator
