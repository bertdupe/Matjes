
module m_H_deriv
!extension of the base type which contains information about how the functional derivatives are to be done
!later also added eval_single, might make to sense to rename this 
use m_H_type, only: t_H_base
use m_deriv_combine
use m_derived_types, only : lattice,number_different_order_parameters
use eval_single, only: t_eval_single
implicit none

private
public t_H
type,abstract,extends(t_H_base) :: t_H
    type(t_deriv_base)          :: deriv(number_different_order_parameters)
    type(t_eval_single)         :: eval_single(number_different_order_parameters)
    logical                     :: deriv_set=.false.
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

    if(this%deriv_set)then
        do i=1,number_different_order_parameters
            Call this%deriv(i)%copy(Hout%deriv(i))
            Call this%eval_single(i)%copy(Hout%eval_single(i))
        enddo
        Hout%deriv_set=.true.
    endif
end subroutine


subroutine set_deriv_single(Ham)
    class(t_H),intent(inout)    :: Ham
    integer                     :: i_mode

    integer         :: Nl,Nr !number of order parameters at left/right side of Hamiltonian 
    integer         :: Nl_order,Nr_order !numbern of order parameters occurances at left/right side of Hamiltonian

    Nl=size(Ham%op_l); Nr=size(Ham%op_r)
    do i_mode=1,number_different_order_parameters
        if(allocated(Ham%deriv(i_mode)%l)) deallocate(Ham%deriv(i_mode)%l)  !shouldn't really be necessary as one could check that is only is called once at the correct moments, but this is easier
        if(allocated(Ham%deriv(i_mode)%r)) deallocate(Ham%deriv(i_mode)%r)
        Nl_order=count(Ham%op_l==i_mode); Nr_order=count(Ham%op_r==i_mode)

        if(Nl+Nr>2)then
            if(Nl_order>0)then
                allocate(Ham%deriv(i_mode)%l,source=t_deriv_l_N(i_mode))
            else
                allocate(t_deriv_l_null::Ham%deriv(i_mode)%l)
            endif
            if(Nr_order>0)then
                allocate(Ham%deriv(i_mode)%r,source=t_deriv_r_N(i_mode))
            else
                allocate(t_deriv_r_null::Ham%deriv(i_mode)%r)
            endif
        elseif(Nl==1.and.Nr==1)then
            if(Nl_order==1.and.Nr_order==1)then
                allocate(t_deriv_l_1_sym::Ham%deriv(i_mode)%l)
                allocate(t_deriv_r_null ::Ham%deriv(i_mode)%r)
            elseif(Nl_order==1)then
                allocate(t_deriv_l_1    ::Ham%deriv(i_mode)%l)
                allocate(t_deriv_r_null ::Ham%deriv(i_mode)%r)
            elseif(Nr_order==1)then
                allocate(t_deriv_l_null ::Ham%deriv(i_mode)%l)
                allocate(t_deriv_r_1    ::Ham%deriv(i_mode)%r)
            else
                allocate(t_deriv_l_null ::Ham%deriv(i_mode)%l)
                allocate(t_deriv_r_null ::Ham%deriv(i_mode)%r)
            endif
       else
            ERROR STOP "THIS SHOULD NOT BE REACHED, something failed setting up the Hamiltonian"
       endif
    enddo

    Ham%deriv_set=.true.
end subroutine

end module
