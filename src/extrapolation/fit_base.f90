module m_fit_base
use m_function_base
use, intrinsic :: iso_fortran_env, only : output_unit
implicit none

type,abstract :: fit
   real(8)                            :: chi                ! final residual difference
   logical                            :: is_set=.false.
   class(function_base),allocatable   :: Func               ! function to be used

   contains
   procedure, NON_OVERRIDABLE :: init_base             ! initialize the parameters
   procedure, NON_OVERRIDABLE :: delete_base           ! delete the parameters

   procedure(int_solve), deferred        :: solve      ! solve the least square problem
   procedure(int_get), deferred          :: get        ! get all parameters of the fit
   procedure(int_execute), deferred      :: execute    ! solve the least square problem
   procedure(int_init), deferred         :: init       ! init all parameters of the fit
end type


abstract interface

    subroutine int_solve(this,f_vec,f_jac)
    import fit
    class(fit), intent(inout)    :: this
    real(8),    intent(inout)    :: f_vec(:)
    real(8),    intent(in)       :: f_jac(:,:)
    end subroutine

    subroutine int_get(this,iter,stop_cr,chi_start,chi_end)
    import fit
    class(fit), intent(in)    :: this
    integer   , intent(inout) :: iter
    real(8)   , intent(inout) :: stop_cr,chi_start,chi_end
    end subroutine

    subroutine int_execute(this,X,Y,tag)
    import fit
    class(fit), intent(inout) :: this
    real(8)   , intent(in)    :: X(:,:),Y(:),tag
    end subroutine

    subroutine int_init(this)
    import fit
    class(fit), intent(inout) :: this
    end subroutine

end interface



private
public :: fit


contains

subroutine init_base(this)
  class(fit), intent(inout) :: this

  this%is_set=.true.

end subroutine

subroutine delete_base(this)
  class(fit), intent(inout) :: this

  this%is_set=.false.

end subroutine

end module
