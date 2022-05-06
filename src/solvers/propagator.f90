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
    use m_vector, only: normalize
    !do the LLG step from data provided in large contiguous array (reshaping to (3,Nsite*Nmag))
    real(8),intent(in)                  ::  damping
    real(8),intent(in),contiguous       ::  M(:,:),B(:,:)
    real(8),intent(inout),contiguous    ::  Mout(:,:)

    !local data
    real(8)     ::  M_norm(size(M,1),size(M,2))
    real(8)     ::  LLG_int(size(M,1),size(M,2))

    if(size(M,1)/=3.or.size(B,1)/=3.or.size(mout,1)/=3) ERROR STOP "LLG INPUT NEEDS TO BE 3-vector"

    M_norm=M
    Call normalize(M_norm)
    Call cross_new(M_norm,B,LLG_int)
    LLG_int=-B-damping*LLG_int
    Call cross_new(M_norm,LLG_int,Mout)
    Mout=Mout/(1.0+damping*damping)

    !ADD EXTERNAL TORQUES IF NECESSARY
end subroutine

!help routine for LLG_new
pure subroutine cross_new(vec1,vec2,vec_out)
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

end module m_propagator
