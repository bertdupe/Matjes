module m_init_random
use m_random_public
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
    integer         :: i,N_rnd

   ! random number
    class(ranbase),allocatable,target   :: random_numbers

    if(mod(dim_mode,3)/=0) STOP "only works for 3-vector like properties"

    call get_ran_type(random_numbers)
    call random_numbers%init_base(size(state))
    call random_numbers%get_extract_list(state)

    state_3(1:3,1:size(state)/3)=>state
    Call normalize(state_3,1.0d-50)

    nullify(state_3)
end subroutine 
end module 
