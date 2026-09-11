module m_fit_linear
use m_fit_base
use, intrinsic :: iso_fortran_env, only : output_unit
use m_invert
implicit none

type,extends(fit) :: fit_linear
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
public :: fit_linear

contains

subroutine execute(this,X,Y,tag)
   class(fit_linear), intent(inout)    :: this
   real(8),    intent(in)    :: X(:,:),Y(:),tag

   ! internal
   real(8),allocatable     :: Jacobian(:,:),moment(:),Gram(:,:),Matrix(:)
   real(8)                 :: rsquare,syx,covariance
   integer                 :: N_y,i,N_variable

   N_y=size(Y)
   N_variable=this%func%N_params
   rsquare=0.0d0
   syx=0.0d0
   covariance=0.0d0

   allocate(Jacobian(N_variable,N_variable),source=0.0d0)
   allocate(Gram(N_variable,N_variable),source=0.0d0)
   allocate(Moment(N_variable),source=0.0d0)
   allocate(Matrix(N_y),source=0.0d0)

   ! calculate the Gram matrix
   Gram=matmul(X,transpose(X))

   ! calculate the Moment matrix
   do i=1,N_variable
      Moment(i)=sum(Y*X(i,:))
   enddo

   call invert(Gram,Jacobian,N_variable)

   this%func%params=matmul(Jacobian,Moment)

   Matrix=matmul(transpose(X),this%func%params)
   rsquare=dot_product(Y-Matrix,Y-Matrix)

   if (this%is_print) call this%print(rsquare,-1.0d0,tag)

end subroutine

subroutine get(this,iter,stop_cr,chi_start,chi_end)
   class(fit_linear), intent(in)    :: this
   integer   , intent(inout) :: iter
   real(8)   , intent(inout) :: stop_cr,chi_start,chi_end

   continue

end subroutine

subroutine solve(this,f_vec,f_jac)
  class(fit_linear), intent(inout)    :: this
  real(8),    intent(inout)    :: f_vec(:)
  real(8),    intent(in)       :: f_jac(:,:)

  continue

end subroutine

subroutine print(this,r2,syx,tag)
class(fit_linear),intent(in)    :: this
real(8),intent(in)       :: r2,syx
real(8),intent(in)       :: tag

   write(output_unit,'(a,f12.6)') 'result of the fit for iteration', tag
   write(output_unit,'(a,E20.12,a)') 'coefficient of determination (r2)', r2, '  (perfect fit r2=0)'
   write(output_unit,'(a,E20.12,a)') 'standard error', syx, '  (perfect fit S=0 - not computed S=1)'

end subroutine

subroutine init(this)
  class(fit_linear), intent(inout)    :: this

  call this%init_base()

end subroutine

end module
