module m_get_random
use mtprng

interface get_rand_classic
 module procedure get_random_standard
 module procedure get_random_Mag
end interface get_rand_classic

interface get_rand_QM
 module procedure get_random_Ising
end interface get_rand_QM
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! functions that returns a random number between 0 and 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function get_random_standard(state) result(Choice)
implicit none
type (mtprng_state),intent(inout) :: state
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
function get_random_Mag(state,N,length)
implicit none
type (mtprng_state),intent(inout) :: state
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
function get_random_Ising(state) result(ising)
implicit none
type (mtprng_state),intent(inout) :: state
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
