!!!!!!!!!!!
! very easy routine that picks up a spin randomly
!!!!
module m_choose_spin
implicit none

interface choose_spin
  module procedure choose_spin_serial
end interface
contains

subroutine choose_spin_serial(iomp,N_spin)
    use m_get_random, only: get_rand_classic
    integer, intent(out)    :: iomp
    integer, intent(in)     :: N_spin
    ! internal variables
    real(8) :: Choice
    !     Choose a random spin place
    choice=get_rand_classic()
    iomp = 1 + NINT(Choice*Dble(N_spin)-0.5d0)
end subroutine choose_spin_serial

end module m_choose_spin
