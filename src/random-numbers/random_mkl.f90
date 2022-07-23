#ifdef CPP_MKL
include 'mkl_vsl.f90'
#endif



module m_random_mkl
use m_random_base
use m_io_files_utils
use m_io_utils


#ifdef CPP_MKL
use MKL_VSL
use MKL_VSL_TYPE
implicit none

private
public :: ranmkl

type,extends(ranbase) :: ranmkl
    type(VSL_STREAM_STATE)  :: stream                     ! Pointer to the stream state structure
    integer                 :: brng = VSL_BRNG_MT19937    ! mersenne twister prng
    integer                 :: method = VSL_RNG_METHOD_GAUSSIAN_BOXMULLER                    ! Generation method ( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER , )
    integer                 :: seed = -717                ! seed for the random number generator

  contains

    procedure  :: init_seed
    procedure  :: get_extract_list
    procedure  :: get_list
    procedure  :: destroy
    procedure  :: read_option
    procedure  :: rand_get
end type

contains

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

    call close_file('input',io_in)

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

    integer(4)  :: stat

    if (mean.gt.1.0d-8) then
       this%mean=mean
       this%is_set=.true.
    endif

    if (.not.this%is_set) return

    stat=vdrnggaussian(this%method,this%stream,this%N,this%x,this%mean,this%sigma)

end subroutine

subroutine rand_get(this,r)
    class(ranmkl),intent(inout)   :: this
    real(8), intent(out)          :: r

    integer(4)  :: stat
    real(8)     :: resu(1)

    stat=vdrnggaussian(this%method,this%stream,1,resu,0.0d0,1.0d0)

    r=resu(1)

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






#endif

end module
