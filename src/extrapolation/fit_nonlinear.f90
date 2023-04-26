module m_fit_nonlinear
use m_fit_base
use, intrinsic :: iso_fortran_env, only : output_unit
use m_io_files_utils
use m_io_utils
implicit none

type,extends(fit) :: fit_nonlinear
   real(8)            :: eps(6)             ! precision of the fit and stoping criterai
   character(len=30)  :: my_fit_name = "internal linear"
   logical            :: is_print = .true.

   contains
   procedure :: solve
   procedure :: get
   procedure :: execute

   procedure :: print
   procedure :: init
end type


private
public :: fit_nonlinear

contains

subroutine execute(this,X,Y,tag)
   class(fit_nonlinear), intent(inout)    :: this
   real(8),intent(in)    :: X(:,:),Y(:),tag

   real(8), allocatable :: F(:),Jacobian(:,:),params(:),delta_params(:)
   integer              :: N_y,i,N_params

   N_y=size(Y)
   N_params=this%func%N_params
   allocate(F(N_y),source=0.0d0)
   allocate(params(N_params),source=this%func%params)
   allocate(delta_params(N_params),source=0.0d0)
   allocate(Jacobian(N_params,N_y),source=0.0d0)

   delta_params=this%eps(1)
   F=this%func%eval_vector([(X(:,i),i=1,N_y)])







   stop 'toto'

end subroutine

subroutine get(this,iter,stop_cr,chi_start,chi_end)
   class(fit_nonlinear), intent(in)    :: this
   integer   , intent(inout) :: iter
   real(8)   , intent(inout) :: stop_cr,chi_start,chi_end

   continue

end subroutine

subroutine solve(this,f_vec,f_jac)
  class(fit_nonlinear), intent(inout)    :: this
  real(8),    intent(inout)    :: f_vec(:)
  real(8),    intent(in)       :: f_jac(:,:)

  continue

end subroutine

subroutine print(this,r2,syx,tag)
  class(fit_nonlinear),intent(in)    :: this
  real(8),intent(in)       :: r2,syx,tag

   write(output_unit,'(a,I10)') 'result of the fit for iteration', tag
   write(output_unit,'(a,E20.12,a)') 'coefficient of determination (r2)', r2, '  (perfect fit r2=1)'
   write(output_unit,'(a,E20.12,a)') 'standard error', syx, '  (perfect fit S=0)'

end subroutine

subroutine init(this)
  class(fit_nonlinear), intent(inout)    :: this

  integer :: io_param

  this%eps=0.0d0
  io_param=open_file_read('input')
  call get_parameter(io_param,'input','eps',this%eps)
  call close_file('input',io_param)

  call this%init_base()

end subroutine

end module
