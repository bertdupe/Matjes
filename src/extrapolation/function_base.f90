module m_function_base
use m_io_files_utils
use m_convert
use m_io_utils
implicit none

type,abstract :: function_base
   real(8), allocatable    :: params(:)          ! parameters of the function
   integer                 :: N_params=-10       ! number of parameters
   logical :: is_set = .False.
   logical :: print = .False.
   character(len=100)  :: name=''
   contains

   procedure, NON_OVERRIDABLE :: init      ! initialize the parameters
   procedure, NON_OVERRIDABLE :: check          ! check all the parameters
   procedure, NON_OVERRIDABLE :: destroy        ! check all the parameters
   procedure, NON_OVERRIDABLE :: print_coef     ! print all the parameters (function dependent)
   procedure, NON_OVERRIDABLE :: init_coef      ! initialize all the parameters (function dependent)

   procedure(show_int),deferred :: show           ! show all the parameters in the outfile (function dependent)
   procedure(set_int),deferred  :: set            ! print all the parameters (function dependent)
   procedure(eval_single_int),deferred :: eval_single    ! evaluate the function
   procedure(eval_vec_int),deferred :: eval_vector    ! evaluate the function

end type


abstract interface

    pure function eval_single_int(this,X) result(Y)
    import function_base
    class(function_base), intent(in) :: this
    real(8), intent(in) :: X
    real(8) :: Y
    end function

    pure function eval_vec_int(this,X) result(Y)
    import function_base
    class(function_base), intent(in) :: this
    real(8), intent(in) :: X(:)
    real(8), allocatable, dimension(:) :: Y
    end function

    subroutine show_int(this,tag)
    import function_base
    class(function_base), intent(in) :: this
    real(8), intent(in) :: tag
    end subroutine

    subroutine set_int(this,N,name)
    import function_base
    class(function_base), intent(inout) :: this
    integer, intent(in) :: N
    character(len=*)  :: name
    end subroutine

end interface

private
public :: function_base

contains



subroutine init(this,N)
  class(function_base), intent(inout) :: this
  integer,intent(in)                  :: N

  if (allocated(this%params)) STOP 'parameters of the function already allocated'

  this%N_params=N
  allocate(this%params(N),source=0.0d0)

end subroutine

subroutine print_coef(this,stats,accesss,tag)
  class(function_base), intent(in) :: this
  real(8), intent(in)   :: tag
  character(len=*), intent(in) :: stats,accesss

  character(len=100) :: fname,form
  integer :: io

  fname=convert('func_',trim(this%name),'_',tag,'.dat')
  form=convert('(',this%N_params,'(E20.12E3,2x))')
  io=open_file_write(fname,stats=stats,accesss=accesss)

  write(io,form) this%params

  call close_file(fname,io)

end subroutine

subroutine destroy(this)
  class(function_base), intent(inout) :: this

  deallocate(this%params)
  this%N_params=-10

end subroutine

function check(this) result(OK)
  class(function_base) :: this
  logical       :: OK

  OK=.false.
  if (allocated(this%params).and.(this%N_params.gt.0).and.(trim(this%name).ne.'')) OK=.true.

end function

subroutine init_coef(this,name)
  class(function_base), intent(inout) :: this
  character(len=*), intent(in) :: name

  character(len=100) :: var_name,fname
  integer :: io_param

  fname='input'
  var_name=convert(trim(name),'_function_params')

  io_param=open_file_read(fname)
  call get_parameter(io_param,trim(fname),trim(var_name),this%params)
  call close_file(fname,io_param)

end subroutine

end module m_function_base
