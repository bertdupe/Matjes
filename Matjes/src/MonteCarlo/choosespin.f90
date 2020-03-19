!!!!!!!!!!!
! very easy routine that picks up a spin randomly
!!!!
module m_choose_spin

interface choose_spin
  module procedure choose_spin_serial
end interface

contains

subroutine choose_spin_serial(iomp,N_spin)
use mtprng
implicit none
integer, intent(inout) :: iomp
integer, intent(in) :: N_spin
! internal variables
type(mtprng_state) :: state
real(kind=8) :: Choice

Choice=0.0d0

!     Choose a random spin place

#ifdef CPP_MRG
 choice=mtprng_rand_real1(state)
#else
 CALL RANDOM_NUMBER(Choice)
#endif
iomp = 1 + NINT(Choice*Dble(N_spin)-0.5d0)


end subroutine choose_spin_serial

end module m_choose_spin
