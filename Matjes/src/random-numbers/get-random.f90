module m_get_random
use mtprng

type(mtprng_state) :: state

interface get_rand_classic
 module procedure get_random_standard
 module procedure get_random_Mag
end interface get_rand_classic

interface get_rand_QM
 module procedure get_random_Ising
end interface get_rand_QM

private
public :: get_rand_QM,get_rand_classic,init_mtprng
contains

subroutine init_mtprng(N)
implicit none
integer, intent(in) :: N

call mtprng_init(N, state)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! functions that returns a random number between 0 and 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function get_random_standard() result(Choice)
implicit none
real(kind=8) :: choice

#ifdef CPP_MRG
      Choice=mtprng_rand_real1(state)
#else
      CALL RANDOM_NUMBER(Choice)
#endif

end function get_random_standard

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! functions that returns a random N dimensions vector if norm length
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function get_random_Mag(N,length)
implicit none
real(kind=8), intent(in) :: length
integer, intent(in) :: N
real(kind=8), dimension(N) :: get_random_Mag(N)
!internal variable
real(kind=8) :: choice,vector(N),norm
integer :: i

get_random_Mag=0.0d0
vector=0.0d0
norm=0.0d0

do i=1,N
#ifdef CPP_MRG
      Choice=mtprng_rand_real1(state)
#else
      CALL RANDOM_NUMBER(Choice)
#endif
      vector(i)=Choice
      norm=norm+Choice**2
enddo

norm=sqrt(norm)
get_random_Mag=vector/norm*length

end function get_random_Mag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! functions that returns a random Ising spin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function get_random_Ising() result(ising)
implicit none
real(kind=8) :: ising
real(kind=8) :: choice

#ifdef CPP_MRG
      Choice=mtprng_rand_real1(state)
#else
      CALL RANDOM_NUMBER(Choice)
#endif

if (choice.lt.0.5d0) then
  ising=-1.0d0
else
  ising=1.0d0
endif

end function get_random_Ising

end module m_get_random
