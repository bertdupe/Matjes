module m_propagator
use m_constants, only : hbar
use m_vector, only : cross,norm,cross_NM
implicit none

private
public :: LLG

contains


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
    Call cross_NM(M_norm,B,LLG_int)
    LLG_int=-B-damping*LLG_int
    Call cross_NM(M_norm,LLG_int,Mout)
    Mout=Mout/(1.0+damping*damping)

    !ADD EXTERNAL TORQUES IF NECESSARY
end subroutine

end module m_propagator
