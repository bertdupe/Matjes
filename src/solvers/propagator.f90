module m_propagator
use m_constants, only : hbar
use m_vector, only : cross,norm,cross_NM
use m_basic_types, only : torque
implicit none

private
public :: LLG

contains


subroutine LLG(B,damping,M,Mout,SOT)
    use m_vector, only: normalize
    !do the LLG step from data provided in large contiguous array (reshaping to (3,Nsite*Nmag))
    real(8),intent(in)                  ::  damping
!    real(8),intent(in)                  ::  P_current(:,:),FL_SOT,DL_SOT
    type(torque),intent(in)             ::  SOT
    real(8),intent(in),contiguous       ::  M(:,:),B(:,:)
    real(8),intent(inout),contiguous    ::  Mout(:,:)

    !local data
    real(8)     ::  M_norm(size(M,1),size(M,2))
    real(8)     ::  LLG_int(size(M,1),size(M,2))
    real(8)     ::  SOT_int(size(M,1),size(M,2)),STT_int(size(M,1),size(M,2))
    real(8)     ::  P_current(size(M,1),size(M,2))
    integer     ::  i


    if(size(M,1)/=3.or.size(B,1)/=3.or.size(mout,1)/=3) ERROR STOP "LLG INPUT NEEDS TO BE 3-vector"

    M_norm=M
    Call normalize(M_norm)

    ! normal equation of motion
    Call cross_NM(M_norm,B,LLG_int)
    LLG_int=B+damping*LLG_int
    Call cross_NM(M_norm,LLG_int,Mout)

    Mout=-Mout/(1.0d0+damping*damping)

    !part that gets the Torques
    !Spin Orbit torques (depends on the polarization of the injected current)
!    if (SOT%is_set) then
!       do i=1,size(M,2)
!          P_current(:,i)=SOT%polarization
!       enddo
!       Call cross_NM(M_norm,P_current,LLG_int)
!       LLG_int=-SOT%FL_damp*P_current-SOT%DL_damp*LLG_int
!       Call cross_NM(M_norm,LLG_int,SOT_int)
!       Mout=Mout+SOT_int
!    endif


    ! Spin transfer torque (depends on the gradient of M)


end subroutine

end module m_propagator
