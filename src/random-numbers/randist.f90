include 'mkl_vsl.f90'

module m_randist
#ifdef CPP_MKL
use MKL_VSL
use MKL_VSL_TYPE
#endif

private
public :: randist

 interface randist
    module procedure gaussianran_mkl
    module procedure gaussianran
    module procedure wiener
 end interface

!interface
!    !manually write interface again since these interfaces seem to differ between different mkl versions
!    function vRngGaussian_fort(method,stream,n,resu,kt,sigma) bind(c, name='vRngGaussian')
!      use iso_c_binding
!      integer(C_int), intent(in)  :: method   ! Generation method ( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER , )
!      type(VSLStreamStatePtr), intent(in)     :: stream   ! Pointer to the stream state structure
!      integer(C_int), intent(in)  :: n        ! Number of random values to be generated.
!      real(C_DOUBLE), intent(in)  :: kt       ! Mean value a.
!      real(C_DOUBLE), intent(in)  :: sigma    ! Standard deviation σ.
!      real(C_DOUBLE), intent(out) :: resu(n)     ! result
!    end function
!
!end interface

contains

!!!! the different ways to generate random forces
real(kind=8) function gaussianran(kt)
implicit none
real(kind=8), intent(in) :: kt
real(kind=8) :: Choice1, Choice2

Choice1=10.0d0
Choice2=10.0d0

  !	open(3,file='unif.dat', access = 'append')
do while (Choice1**2+Choice2**2.gt.1.0d0)
  call RANDOM_NUMBER(Choice1)
  call RANDOM_NUMBER(Choice2)
!write(3,*) Choice1, " " ,Choice2
  		
  		
  Choice1=(Choice1*2.0d0-1.0d0)
  Choice2=(Choice2*2.0d0-1.0d0)
enddo

! close(3)

!
! according to the PhD of Schieback (2010) from Constanz p. 35
!

!gaussianran=sqrt(2.0d0*kT)*Choice1*sqrt(-2.0d0*dlog(Choice1**2+Choice2**2)/(Choice1**2+Choice2**2))
gaussianran=sqrt(kT)*Choice1*sqrt(-2.0d0*dlog(Choice1**2+Choice2**2)/(Choice1**2+Choice2**2))
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

#ifdef CPP_MKL
!!!! Gaussian random number generator of MKL
function gaussianran_mkl(kt,sigma)
implicit none
real(8), intent(in)     :: kt,sigma
integer(4)              :: stat
type(VSL_STREAM_STATE)  :: stream
integer                 :: brng,method,seed,n
real(8)                 :: gaussianran_mkl(1)

brng=VSL_BRNG_MT19937
method=VSL_RNG_METHOD_GAUSSIAN_BOXMULLER
seed=7

    ! initialization
    stat=vslnewstream( stream, brng,  seed )

    stat=vdrnggaussian(method,stream,1,gaussianran_mkl,kt,sigma)

    ! destroy
    stat=vsldeletestream( stream )

end function

#endif

end module m_randist
