module m_randist
 interface randist
    module procedure gaussianran
    module procedure wiener
 end interface

contains
real(kind=8) function gaussianran(kt)
use m_dynamic
implicit none
real(kind=8), intent(in) :: kt
real(kind=8) :: Choice1, Choice2

call RANDOM_NUMBER(Choice1)
call RANDOM_NUMBER(Choice2)

Choice1=(Choice1*2.0d0-1.0d0)/dsqrt(2.0d0)
Choice2=(Choice2*2.0d0-1.0d0)/dsqrt(2.0d0)

!      if (dlog(Choice1**2+Choice2**2)/(Choice1**2+Choice2**2).gt.0.0d0) STOP

gaussianran=dsqrt(2.0d0*kT)*Choice1*dsqrt(-dlog(Choice1**2+Choice2**2)/(Choice1**2+Choice2**2))

! as expressed in Mentink et al J. Phys. C 22, 176001 (2010)
!      gaussianran=2.0d0*damping*kT/(1+damping**2)*Choice1*dsqrt(-2.0d0*dlog( &
!         Choice1**2+Choice2**2)/(Choice1**2+Choice2**2))
end function

      real(kind=8) function wiener(state)
      use m_dynamic, only : timestep
      use mtprng
      implicit none
      type(mtprng_state), intent(inout) :: state
!dummy
      real(kind=8) :: Choice

#ifdef CPP_MRG
      Choice=mtprng_rand_real1(state)
#else
      CALL RANDOM_NUMBER(Choice)
#endif

      wiener=sqrt(timestep)*(Choice*2.0d0-1.0d0)

      end function wiener

end module m_randist
