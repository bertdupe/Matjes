
module m_H_deriv
!extension of the base type which contains information about how the functional derivatives are to be done
!later also added eval_single, might make to sense to rename this 
use m_H_type, only: t_H_base
use m_deriv_public, only: t_deriv, set_deriv
use m_derived_types, only : lattice,number_different_order_parameters
use eval_single, only: t_eval_single
implicit none

private
public t_H
type,abstract,extends(t_H_base) :: t_H
    type(t_deriv)               :: deriv(number_different_order_parameters)
    type(t_eval_single)         :: eval_single(number_different_order_parameters)
contains
    procedure   :: copy_deriv
    procedure   :: set_deriv => set_deriv_single
    procedure   :: finish_setup
end type

contains


subroutine finish_setup(this)
    class(t_H),intent(inout)    :: this

    integer     :: i_mode

    Call this%finish_setup_base()
    Call this%set_deriv()
    do i_mode=1,number_different_order_parameters
        Call this%eval_single(i_mode)%set(this,i_mode)
    enddo
end subroutine

subroutine copy_deriv(this,Hout)
    class(t_H),intent(in)     :: this
    class(t_H),intent(inout)  :: Hout

    integer ::  i

    do i=1,number_different_order_parameters
        Call this%deriv(i)%copy(Hout%deriv(i))
        Call this%eval_single(i)%copy(Hout%eval_single(i))
    enddo
end subroutine


subroutine set_deriv_single(Ham)
    class(t_H),intent(inout)    :: Ham
    integer                     :: i
    do i=1,number_different_order_parameters
        Call set_deriv(Ham%deriv(i),i,Ham%op_l,Ham%op_r)      
    enddo
end subroutine

end module
