module m_H_deriv
use m_H_type, only: t_H_base
use m_deriv_combine
use m_derived_types, only : lattice,number_different_order_parameters

private
public t_H
type,abstract,extends(t_H_base) :: t_H
    type(t_deriv_base)          :: deriv(number_different_order_parameters)
    logical                     :: deriv_set=.false.
contains
    procedure   :: copy_deriv
    procedure   :: set_deriv => set_deriv_single
end type

contains

subroutine copy_deriv(this,Hout)
    class(t_H),intent(in)     :: this
    class(t_H),intent(inout)  :: Hout

    integer ::  i

    if(this%deriv_set)then
        do i=1,number_different_order_parameters
            Call this%deriv(i)%copy(Hout%deriv(i))
        enddo
        Hout%deriv_set=.true.
    endif
end subroutine


subroutine set_deriv_single(Ham)
    class(t_H),intent(inout)    :: Ham
    integer                     :: i_mode

    integer         :: Nl,Nr !number of order parameters at left/right side of Hamiltonian 
    integer         :: Nl_order,Nr_order !numbern of order parameters occurances at left/right side of Hamiltonian

    do i_mode=1,number_different_order_parameters
        Nl=size(Ham%op_l); Nr=size(Ham%op_r)
        Nl_order=count(Ham%op_l==i_mode); Nr_order=count(Ham%op_r==i_mode)

        if(Nl==1.and.Nr==1.and.Nl_order==1.and.Nr_order==1)then
            allocate(t_deriv_l_1_sym::Ham%deriv(i_mode)%l)
            allocate(t_deriv_r_null::Ham%deriv(i_mode)%r)
            cycle
        endif

        if(Nl_order>0)then
            if(Nl==1)then
                allocate(t_deriv_l_1::Ham%deriv(i_mode)%l)
            elseif(Nl>1)then
                allocate(Ham%deriv(i_mode)%l,source=t_deriv_l_N(i_mode))
            endif
        else
            allocate(t_deriv_l_null::Ham%deriv(i_mode)%l)
        endif

        if(Nr_order>0)then
            if(Nr==1)then
                allocate(t_deriv_r_1::Ham%deriv(i_mode)%r)
            elseif(Nr>1)then
                allocate(Ham%deriv(i_mode)%r,source=t_deriv_r_N(i_mode))
            endif
        else
            allocate(t_deriv_r_null::Ham%deriv(i_mode)%r)
        endif
    enddo
    Ham%deriv_set=.true.
end subroutine

end module
