module m_randist

 interface randist
    module procedure gaussianran
    module procedure wiener
 end interface

private
public :: randist

contains

!!!! the different ways to generate random forces
real(kind=8) function gaussianran(kt)
implicit none
real(kind=8), intent(in) :: kt
real(kind=8) :: Choice1, Choice2

Choice1=10.0d0
Choice2=10.0d0

do while (Choice1**2+Choice2**2.gt.1.0d0)
  call RANDOM_NUMBER(Choice1)
  call RANDOM_NUMBER(Choice2)
  Choice1=(Choice1*2.0d0-1.0d0)
  Choice2=(Choice2*2.0d0-1.0d0)
enddo

!
! according to the PhD of Schieback (2010) from Constanz p. 35
!

gaussianran=sqrt(2.0d0*kT)*Choice1*sqrt(-2.0d0*dlog(Choice1**2+Choice2**2)/(Choice1**2+Choice2**2))

end function


!!!! wiener process
real(kind=8) function wiener()
use mtprng
implicit none
!dummy
type(mtprng_state) :: state
real(kind=8) :: Choice

#ifdef CPP_MRG
Choice=mtprng_rand_real1(state)
#else
CALL RANDOM_NUMBER(Choice)
#endif

wiener=(Choice*2.0d0-1.0d0)

end function wiener

end module m_randist
