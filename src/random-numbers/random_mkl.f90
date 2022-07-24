#ifdef CPP_MKL
include 'mkl_vsl.f90'
#endif



module m_random_mkl
use m_random_base
use m_io_files_utils
use m_io_utils
use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit, output_unit


#ifdef CPP_MKL
use MKL_VSL
use MKL_VSL_TYPE
implicit none

private
public :: ranmkl

type,extends(ranbase) :: ranmkl
    type(VSL_STREAM_STATE)  :: stream                     ! Pointer to the stream state structure
    integer                 :: brng = VSL_BRNG_MT19937    ! mersenne twister prng
    integer                 :: method = VSL_RNG_METHOD_UNIFORM_STD_ACCURATE                    ! Generation method
    integer                 :: seed = -717                ! seed for the random number generator

    procedure(int_rnd_distrib),pointer,pass,public :: rnd_distrib => null()

  contains

    procedure  :: init_seed
    procedure  :: get_extract_list
    procedure  :: get_list
    procedure  :: destroy
    procedure  :: read_option
    procedure  :: rand_get

end type

abstract interface
    subroutine int_rnd_distrib(this)
       import :: ranmkl
       class(ranmkl), intent(inout)     :: this
    end subroutine
end interface


contains

subroutine read_option(this)
    class(ranmkl),intent(inout)   :: this

    integer :: io_in

    io_in=open_file_read('input')

    call get_parameter(io_in,'input','vsl_method',this%method)
    call get_parameter(io_in,'input','vsl_brng',this%brng)
    call get_parameter(io_in,'input','print_rnd',this%print)
    call get_parameter(io_in,'input','rnd_noise',this%is_set)
    call get_parameter(io_in,'input','mean',this%mean)
    call get_parameter(io_in,'input','sigma',this%sigma)
    call get_parameter(io_in,'input','min_rnd_val',this%min_rnd_val)
    call get_parameter(io_in,'input','max_rnd_val',this%max_rnd_val)
    call get_parameter(io_in,'input','name_rnd',this%name)
    call close_file('input',io_in)

    call select_rnd_in_lib(this%name,this%rnd_distrib,this%method)

end subroutine

!!!! initialize the random number generator of MKL
subroutine init_seed(this)
    class(ranmkl),intent(inout)   :: this

    integer(4)                    :: stat

    ! initialization
    stat=vslnewstream( this%stream, this%brng,  this%seed )

end subroutine

!!!! destroy the random number generator of MKL
subroutine destroy(this)
    class(ranmkl),intent(in)   :: this

    integer(4)                 :: stat

    ! destroy
    stat=vsldeletestream( this%stream )

end subroutine

!!!! Gaussian random number generator of MKL
subroutine get_list(this,mean)
    class(ranmkl),intent(inout)  :: this
    real(8),intent(in)           :: mean

    if (mean.gt.1.0d-8) then
       this%mean=mean
       this%is_set=.true.
    endif

    if (.not.this%is_set) return

    call this%rnd_distrib()

end subroutine

subroutine rand_get(this,r)
    class(ranmkl),intent(inout)   :: this
    real(8), intent(out)          :: r

    integer(4)  :: stat
    real(8)     :: res(1)

    stat=vdRngUniform(this%method,this%stream,1,res,this%min_rnd_val,this%max_rnd_val)

    r=res(1)

end subroutine

!!!! Gaussian random number generator of MKL
subroutine get_extract_list(this,resu)
    class(ranmkl),intent(inout)         :: this
    real(8), intent(inout),allocatable  :: resu(:)

    if (size(resu).ne.this%N) stop 'size of the result does not match the number of randome numbers'
    if (allocated(resu)) stop 'resu matrix already allocated'
    allocate(resu(this%N),source=0.0d0)

    call this%get_list(this%mean)

    call this%extract_list(resu)

end subroutine



















subroutine select_rnd_in_lib(name,rnd_distrib,method)
character(len=3),intent(in)   :: name
integer,intent(inout)         :: method
procedure(int_rnd_distrib),pointer,intent(out)  :: rnd_distrib

select case (name)
   case('uni')
      write(output_unit,'(/a/)') "uniform distribution selected"
      rnd_distrib => rand_uniform
      method=VSL_RNG_METHOD_UNIFORM_STD_ACCURATE

   case('nor')
      write(output_unit,'(/a)') "normal (Gaussian) distribution"
      write(output_unit,'(a/)') "mean => displacement, sigma => Standard deviation"
      rnd_distrib => rand_normal
      method=VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2

   case('exp')
      write(output_unit,'(/a)') "exponential distribution selected"
      write(output_unit,'(a/)') "mean => displacement, sigma => Scalefactor"
      rnd_distrib => rand_exponential
      method=VSL_RNG_METHOD_EXPONENTIAL_ICDF_ACCURATE

   case('gam')
      write(output_unit,'(/a)') "gamma distribution selected"
      write(output_unit,'(a/)') "min_rnd_val => Shape, mean => Displacement ,sigma => Scalefactor"
      rnd_distrib => rand_gamma
      method=VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE

   case('chi')
      write(output_unit,'(/a)') "chi square distribution selected"
      write(output_unit,'(a/)') "mean => degrees of freedom"
      rnd_distrib => rand_chi_square
      method=VSL_RNG_METHOD_CHISQUARE_CHI2GAMMA

   case('wei')
      write(output_unit,'(/a)') "Weibull distribution selected"
      write(output_unit,'(a/)') "mean => Shape, min_rnd_val => displacement, sigma => Scalefactor"
      rnd_distrib => rand_weibull
      method=VSL_RNG_METHOD_WEIBULL_ICDF_ACCURATE

   case('cau')
      write(output_unit,'(/a)') "Cauchy distribution selected"
      write(output_unit,'(a/)') "mean => Displacement, sigma => Scalefactor"
      rnd_distrib => rand_cauchy
      method=VSL_RNG_METHOD_CAUCHY_ICDF

   case('lap')
      write(output_unit,'(/a)') "laplace distribution selected"
      write(output_unit,'(a/)') "mean => mean value, sigma => Scalefactor"
      rnd_distrib => rand_laplace
      method=VSL_RNG_METHOD_LAPLACE_ICDF

   case('log')
      write(output_unit,'(/a)') "log-normal distribution selected"
      write(output_unit,'(a/)') "mean => Average, sigma=> standard deviation, min_rnd_val=> Displacement, max_rnd_val => Scalefactor"
      rnd_distrib => rand_log_normal
      method=VSL_RNG_METHOD_LOGNORMAL_BOXMULLER2_ACCURATE

   case('bet')
      write(output_unit,'(/a)') "Beta distribution selected"
      write(output_unit,'(a/)') "min_rnd_val => Shape p, max_rnd_val => Shape q, mean => displacement, sigma => scale factor"
      rnd_distrib => rand_beta
      method=VSL_RNG_METHOD_BETA_CJA_ACCURATE

   case('ray')
      write(output_unit,'(/a)') "Rayleigh distribution selected"
      write(output_unit,'(a/)') "mean => Displacement, sigma => Scalefactor"
      rnd_distrib => rand_rayleigh
      method=VSL_RNG_METHOD_RAYLEIGH_ICDF_ACCURATE

   case('gum')
      write(output_unit,'(/a)') "Gumbel distribution selected"
      write(output_unit,'(a/)') "mean => Displacement, sigma => Scalefactor"
      rnd_distrib => rand_gumbel
      method=VSL_RNG_METHOD_GUMBEL_ICDF

       case default
         stop 'please select a distribution for your random numbers'
    end select

end subroutine










!
! uniform random Number Generators in Fortran
!
subroutine rand_uniform(this)
class(ranmkl),intent(inout)       :: this
integer(4)                        :: stat

   stat=vdRngUniform(this%method,this%stream,this%N,this%x,this%min_rnd_val,this%max_rnd_val)

end subroutine

!
! Random Sample from normal (Gaussian) distribution
!
!interface
!    !manually write interface again since these interfaces seem to differ between different mkl versions
!    function vRngGaussian_fort(method,stream,n,resu,mean,sigma) bind(c, name='vRngGaussian')
!      use iso_c_binding
!      integer(C_int), intent(in)  :: method   ! Generation method ( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER , )
!      type(VSLStreamStatePtr), intent(in)     :: stream   ! Pointer to the stream state structure
!      integer(C_int), intent(in)  :: n        ! Number of random values to be generated.
!      real(C_DOUBLE), intent(in)  :: a        ! Mean value
!      real(C_DOUBLE), intent(in)  :: b        ! Standard deviation σ.
!      real(C_DOUBLE), intent(out) :: resu(n)     ! result
!    end function
!
!end interface

subroutine rand_normal(this)
class(ranmkl),intent(inout) :: this
integer(4)                   :: stat

   stat=vdRngGaussian(this%method,this%stream,this%N,this%x,this%mean,this%sigma)

end subroutine

!
! Random sample from an exponential distribution
!
subroutine rand_exponential(this)
class(ranmkl),intent(inout) :: this
integer(4)                   :: stat

   stat=vdRngExponential(this%method,this%stream,this%N,this%x,this%mean,this%sigma)

end subroutine

!
! Return a random sample from a gamma distribution
!
subroutine rand_gamma(this)
class(ranmkl),intent(inout) :: this
integer(4)                   :: stat

   stat=vdRngGamma(this%method,this%stream,this%N,this%x,this%min_rnd_val,this%mean,this%sigma)

end subroutine

!
! ## return a random sample from a chi square distribution
! ## with the specified degrees of freedom
!
subroutine rand_chi_square(this)
class(ranmkl),intent(inout) :: this
integer(4)                   :: stat

   stat=vdRngChiSquare(this%method,this%stream,this%N,this%x,int(this%mean))

end subroutine

!
!## return a sample from a Weibull distribution
!
subroutine rand_weibull(this)
class(ranmkl),intent(inout) :: this
integer(4)                   :: stat

   stat=vdRngWeibull(this%method,this%stream,this%N,this%x,this%mean,this%min_rnd_val,this%sigma)

end subroutine

!
!## return a random sample from a Cauchy distribution
!
subroutine rand_cauchy(this)
class(ranmkl),intent(inout) :: this
integer(4)                   :: stat

   stat=vdRngCauchy(this%method,this%stream,this%N,this%x,this%mean,this%sigma)

end subroutine

!
!## return a random sample from a Laplace distribution
!## The Laplace distribution is also known as the double exponential distribution.
!
subroutine rand_laplace(this)
class(ranmkl),intent(inout) :: this
integer(4)                   :: stat

   stat=vdRngLaplace(this%method,this%stream,this%N,this%x,this%mean,this%sigma)

end subroutine

!
! ## return a random sample from a log-normal distribution
!
subroutine rand_log_normal(this)
class(ranmkl),intent(inout) :: this
integer(4)                   :: stat

   stat=vdRngLognormal(this%method,this%stream,this%N,this%x,this%mean,this%sigma,this%min_rnd_val,this%max_rnd_val)

end subroutine

!
! ## return a random sample from a beta distribution
!
subroutine rand_beta(this)
class(ranmkl),intent(inout) :: this
integer(4)                   :: stat

   stat=vdRngBeta(this%method,this%stream,this%N,this%x,this%min_rnd_val,this%max_rnd_val,this%mean,this%sigma)

end subroutine

subroutine rand_rayleigh(this)
class(ranmkl),intent(inout) :: this
integer(4)                   :: stat

   stat=vdRngRayleigh(this%method,this%stream,this%N,this%x,this%mean,this%sigma)

end subroutine

subroutine rand_gumbel(this)
class(ranmkl),intent(inout) :: this
integer(4)                   :: stat

   stat=vdRngGumbel(this%method,this%stream,this%N,this%x,this%mean,this%sigma)

end subroutine



#endif

end module
