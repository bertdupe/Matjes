module m_random_prng
use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit, output_unit, int64
use m_random_base
use mtprng
use m_io_files_utils
use m_io_utils
use m_random_number_library
implicit none

private
public :: ranprng

type,extends(ranbase) :: ranprng

  procedure(int_rnd_distrib),pointer,nopass  :: rnd_distrib => null()
  type(mtprng_state)  ::  state

  contains

    procedure  :: init_seed
    procedure  :: get_extract_list
    procedure  :: get_list
    procedure  :: destroy
    procedure  :: read_option
    procedure  :: rand_get
end type

contains

subroutine read_option(this)
    class(ranprng), intent(inout)  :: this

    integer :: io_in

    io_in=open_file_read('input')
    call get_parameter(io_in,'input','print_rnd',this%print)
    call get_parameter(io_in,'input','rnd_noise',this%is_set)
    call get_parameter(io_in,'input','mean',this%mean)
    call get_parameter(io_in,'input','sigma',this%sigma)
    call get_parameter(io_in,'input','name_rnd',this%name)
    call get_parameter(io_in,'input','min_rnd_val',this%min_rnd_val)
    call get_parameter(io_in,'input','max_rnd_val',this%max_rnd_val)
    call close_file('input',io_in)

    call select_rnd_in_lib(this%name,this%rnd_distrib)

end subroutine

!!!! initialize the random number generator of GNU
subroutine init_seed(this)
    class(ranprng),intent(inout)   :: this

   integer, allocatable :: seed(:)
   integer :: i, n, un, istat, dt(8), pid, getpid
   integer(int64) :: t

   call random_seed(size = n)
   allocate(seed(n))
   ! First try if the OS provides a random number generator
   open(newunit=un, file="/dev/urandom", access="stream", &
    form="unformatted", action="read", status="old", iostat=istat)
   if (istat == 0) then
      read(un) seed
      close(un)
   else
   ! Fallback to XOR:ing the current time and pid. The PID is
   ! useful in case one launches multiple instances of the same
   ! program in parallel.
      call system_clock(t)
      if (t == 0) then
         call date_and_time(values=dt)
         t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24_int64 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
      end if
      pid = getpid()
      t = ieor(t, int(pid, kind(t)))
      do i = 1, n
         seed(i) = lcg(t)
      end do
   end if
   call mtprng_init_by_array(seed,this%state)

end subroutine

!!!! destroy the random number generator of GNU
subroutine destroy(this)
    class(ranprng),intent(in)   :: this

    continue

end subroutine

subroutine get_list(this,mean)
    class(ranprng),intent(inout)  :: this
    real(8),intent(in)            :: mean

    integer :: i
    real(8) :: res

    if (mean.gt.1.0d-8) then
       this%mean=mean
       this%is_set=.true.
    endif

    if (.not.this%is_set) return

    do i=1,this%N
       res=this%rnd_distrib(this)
       this%x(i)=res
    enddo

end subroutine


subroutine rand_get(this,r)
    class(ranprng),intent(inout)  :: this
    real(8), intent(out)          :: r

    r=mtprng_rand_real1(this%state)

end subroutine

subroutine get_extract_list(this,resu)
    class(ranprng),intent(inout)        :: this
    real(8), intent(inout),allocatable  :: resu(:)

    if (size(resu).ne.this%N) stop 'size of the result does not match the number of randome numbers'
    if (allocated(resu)) stop 'resu matrix already allocated'
    allocate(resu(this%N),source=0.0d0)

    call this%get_list(this%mean)

    call this%extract_list(resu)

end subroutine




! This simple PRNG might not be good enough for real work, but is
! sufficient for seeding a better PRNG.
function lcg(s)
integer :: lcg
integer(int64) :: s
if (s == 0) then
   s = 104729
else
   s = mod(s, 4294967296_int64)
end if
s = mod(s * 279470273_int64, 4294967291_int64)
lcg = int(mod(s, int(huge(0), int64)), kind(0))
end function lcg

end module
