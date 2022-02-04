module m_randist
use mtprng
use m_random_base
use iso_fortran_env, only: int64
use m_io_files_utils
use m_io_utils
implicit none

private
public :: ranint,randist

type,extends(ranbase) :: ranint

  contains

    procedure  :: init_seed
    procedure  :: get_list
    procedure  :: destroy
    procedure  :: get_extract_list
    procedure  :: read_option

end type

 interface randist
    module procedure gaussianran
    module procedure wiener
 end interface

contains


subroutine get_list(this,a,b)
class(ranint), intent(inout)  :: this
real(8), intent(in)           :: a,b

integer :: i

do i=1,this%N
   this%x(i)=randist(a)
enddo

end subroutine

subroutine get_extract_list(this,a,b,resu)
class(ranint), intent(inout)  :: this
real(8), intent(in)           :: a,b
real(8), intent(inout)        :: resu(:)

resu=0.0d0

call this%get_list(a,b)

call this%extract_list(resu)

end subroutine

subroutine destroy(this)
class(ranint), intent(in)  :: this

return

end subroutine


subroutine init_seed(this)
class(ranint), intent(in)  :: this

integer, allocatable :: seed(:)
integer :: i, n, un, istat, dt(8), pid
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
   call random_seed(put=seed)

end subroutine init_seed

!!!! the different ways to generate random forces
real(kind=8) function gaussianran(kt)
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





subroutine read_option(this)
    class(ranint), intent(inout)  :: this

    integer :: io_in

    io_in=open_file_read('input')
    call get_parameter(io_in,'input','print_rnd',this%print)
    call close_file('input',io_in)

end subroutine

end module m_randist
