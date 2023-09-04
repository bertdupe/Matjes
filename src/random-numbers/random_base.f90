module m_random_base
use m_io_files_utils
use m_convert
use, intrinsic :: iso_fortran_env, only : error_unit,output_unit
use mpi_basic
#ifdef CPP_MPI
use mpi
#endif
implicit none

private
public :: ranbase

type, abstract :: ranbase
    real(8),allocatable          :: x(:)      ! variable containing all the random numbers
    integer                      :: N         ! number of random numbers to be generated
    logical                      :: print = .false.   ! print out the numbers
    logical                      :: io_read = .false.   ! print out the numbers
    character(len=30)            :: name='uniform'
    real(8)                      :: mean=0.0d0        ! random numbers are generated with a mean value (by default 0,5 so mean*[0,1])
    real(8)                      :: sigma=1.0d0       ! spread of generated numbers
    real(8)                      :: max_rnd_val=1.0d0
    real(8)                      :: min_rnd_val=0.0d0
    logical                      :: is_set = .false.

contains
    ! defered type
    procedure(int_seed),deferred                   :: init_seed
    procedure(int_getx),deferred                   :: get_list
    procedure(int_destroy),deferred                :: destroy
    procedure(int_gextra_list),deferred            :: get_extract_list
    procedure(int_rw_option),deferred              :: read_option
    procedure(int_rand_get),deferred               :: rand_get

    ! base function
    procedure, NON_OVERRIDABLE     :: init_base
    procedure, NON_OVERRIDABLE     :: print_base
    procedure, NON_OVERRIDABLE     :: read_base
    procedure, NON_OVERRIDABLE     :: extract_list
    procedure, NON_OVERRIDABLE     :: extract_size
    procedure, NON_OVERRIDABLE     :: bcast_val
    procedure, NON_OVERRIDABLE     :: set
end type


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
abstract interface

    subroutine int_rand_get(this,r)
       import ranbase
       class(ranbase), intent(inout)   :: this
       real(kind=8)  , intent(out)     :: r
    end subroutine

    subroutine int_seed(this)
       import ranbase
       class(ranbase), intent(inout)     :: this
    end subroutine

    subroutine int_getx(this,kt)
       import ranbase
       class(ranbase), intent(inout)  :: this
       real(8), intent(in)            :: kt
    end subroutine

    subroutine int_destroy(this)
       import ranbase
       class(ranbase), intent(in)     :: this
    end subroutine

    subroutine int_gextra_list(this,resu)
       import ranbase
       class(ranbase), intent(inout)      :: this
       real(8),intent(inout)              :: resu(:)
    end subroutine

    subroutine int_rw_option(this)
       import ranbase
       class(ranbase), intent(inout)  :: this
    end subroutine

end interface


contains

subroutine set(this,mean,sigma,max_rnd_val,min_rnd_val)
    class(ranbase),intent(inout) :: this
    real(8), intent(in),optional :: mean,sigma,max_rnd_val,min_rnd_val

    logical :: changed

    changed=.false.

    if (present(mean)) then
       if (abs(this%mean-mean).gt.epsilon(0.0d0)) then
          this%mean=mean
          changed=.true.
       endif
    endif

    if (present(sigma)) then
       if (abs(this%sigma-sigma).gt.epsilon(0.0d0)) then
          this%sigma=sigma
          changed=.true.
       endif
    endif

    if (present(max_rnd_val)) then
       if (abs(this%max_rnd_val-max_rnd_val).gt.epsilon(0.0d0)) then
          this%max_rnd_val=max_rnd_val
       endif
    endif

    if (present(min_rnd_val)) then
       if (abs(this%min_rnd_val-min_rnd_val).gt.epsilon(0.0d0)) then
          this%min_rnd_val=min_rnd_val
       endif
    endif

    if (changed) then
       write(output_unit,'(/a,2x,10a)') "parameter for random number generator:",this%name
       write(output_unit,'(4(a,x,f10.8,x)/)') "mean", this%mean, "sigma", this%sigma, "max_rnd_val", this%max_rnd_val, "min_rnd_val", this%min_rnd_val
    endif

end subroutine

subroutine extract_list(this,x)
    class(ranbase),intent(in)    :: this
    real(8), intent(inout)       :: x(:)

    integer :: N

    N=size(x)
    x=this%x(1:N)

end subroutine

subroutine extract_size(this,N)
    class(ranbase),intent(in)    :: this
    integer, intent(out)         :: N

    N=this%N

end subroutine

subroutine init_base(this,N)
    class(ranbase),intent(inout) :: this
    integer, intent(in)          :: N

    if (allocated(this%x)) STOP "cannot allocate random number vector - already allocated"

    this%N=N
    allocate(this%x(N),source=0.0d0)

end subroutine


subroutine print_base(this,tag)
    class(ranbase),intent(in) :: this
    integer, intent(in)          :: tag

    integer :: i,io_out
    character(len=100) :: fname

    if (.not.allocated(this%x)) STOP "cannot print random numbers vector - list not allocated"

    fname=convert('rnd_num_',tag)
    io_out=open_file_write(trim(fname))
    do i=1,this%N
       write(io_out,'(E20.12E3)') this%x(i)
    enddo
    call close_file(fname,io_out)

end subroutine

subroutine read_base(this)
    class(ranbase),intent(inout) :: this

    integer :: i,io_out
    character(len=100) :: fname

    if (.not.allocated(this%x)) STOP "cannot read random numbers vector - list not allocated"

    fname='rnd_num_in'
    io_out=open_file_read(trim(fname))
    do i=1,this%N
       write(io_out,'(E20.12E3)') this%x(i)
    enddo
    call close_file(fname,io_out)

end subroutine

subroutine bcast_val(this,comm)
    class(ranbase),intent(inout)    ::  this
    type(mpi_type),intent(in)       ::  comm

    integer     ::  ierr

#ifdef CPP_MPI

Call MPI_Bcast(this%x , this%N , MPI_REAL8 , comm%mas , comm%com ,ierr)

#else
    continue
#endif
end subroutine

end module
