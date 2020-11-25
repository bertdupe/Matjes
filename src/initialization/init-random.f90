module m_init_random
implicit none
private
public :: init_random

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the starting configuration as a random spin configuration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_random(dim_mode,state)
    use m_vector,only: normalize
    integer, intent(in)             :: dim_mode
    real(8),pointer,intent(inout)   :: state(:)

    real(8),pointer :: state_3(:,:)
    integer         :: i

    if(mod(dim_mode,3)/=0) STOP "only works for 3-vector like properties"
#ifdef CPP_MRG
    STOP "cannot use init_random with CPP_MRG, needs to be implemented"
#else
    CALL RANDOM_NUMBER(state)
#endif
    
    state_3(1:3,1:size(state)/3)=>state
    Call normalize(state_3,1.0d-50)
    !do i=1,size(state_3,2)
    !    state_3(:,i)=(get_rand_classic(3,1.0d0)-0.5d0)*2.0d0
    !enddo
    nullify(state_3)
end subroutine 
end module 
