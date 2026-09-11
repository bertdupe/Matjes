module m_func_linear
use m_function_base
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit
implicit none

type,extends(function_base) :: func_linear

   contains
   procedure :: show    ! show all the parameters in the outfile (function dependent)
   procedure :: set     ! print all the parameters (function dependent)
   procedure :: eval_single
   procedure :: eval_vector
end type

private
public :: func_linear

contains

pure function eval_single(this,X) result(Y)
  class(func_linear), intent(in) :: this
  real(8), intent(in) :: X
  real(8)             :: Y

  real(8) :: a,b

  a=this%params(1)
  b=this%params(2)

  Y=a*X+b

end function

pure function eval_vector(this,X) result(Y)
  class(func_linear), intent(in) :: this
  real(8), intent(in) :: X(:)
  real(8), allocatable, dimension(:) :: Y

  real(8) :: a,b
  integer :: N

  a=this%params(1)
  b=this%params(2)
  N=size(X)
  allocate(Y(N),source=0.0d0)

  Y=a*X+b

end function

subroutine set(this,N,name)
  class(func_linear), intent(inout) :: this
  integer, intent(in) :: N
  character(len=*)  :: name

  call this%init(N)
  this%name=name
  call this%init_coef(name)

  if (this%check()) then
     write(output_unit,'(a)') 'linear function is set'
     this%is_set=.true.
  else
     stop 'linear function is NOT set properly'
  endif

end subroutine

subroutine show(this,tag)
  class(func_linear), intent(in) :: this
  real(8), intent(in) :: tag

  write(output_unit,'(2a,a,f12.6)') 'parameters for   ', this%name, 'step:   ', tag
  write(output_unit,'(a,2(a,f12.6))') 'ax+b  ','a=',this%params(1),'  b=',this%params(2)

end subroutine

end module
