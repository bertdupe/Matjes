module m_deriv_public
!use m_derived_types,only : lattice, number_different_order_parameters
use m_H_deriv, only: t_H
!use m_deriv_base, only: t_deriv_base
!
!use m_deriv_null
!use m_deriv_rank2
!use m_deriv_rankN


private
public :: set_deriv!, t_deriv_base
contains


!subroutine set_deriv_single(Ham)
!    class(t_H) :: Ham
!    integer         :: i_mode
!
!    integer         :: Nl,Nr !number of order parameters at left/right side of Hamiltonian 
!    integer         :: Nl_order,Nr_order !numbern of order parameters occurances at left/right side of Hamiltonian
!
!    do i_mode=1,number_different_order_parameters
!        Nl=size(Ham%op_l); Nr=size(Ham%op_r)
!        Nl_order=count(Ham%op_l==i_mode); Nr_order=count(Ham%op_r==i_mode)
!
!        if(Nl==1.and.Nr==1.and.Nl_order==1.and.Nr_order==1)then
!            allocate(t_deriv_l_1_sym::Ham%deriv(i_mode)%l)
!            allocate(t_deriv_r_null::Ham%deriv(i_mode)%r)
!            cycle
!        endif
!
!        if(Nl_order>0)then
!            if(Nl==1)then
!                allocate(t_deriv_l_1::Ham%deriv(i_mode)%l)
!            elseif(Nl>1)then
!                allocate(Ham%deriv(i_mode)%l,source=t_deriv_l_N(i_mode))
!            endif
!        else
!            allocate(t_deriv_l_null::Ham%deriv(i_mode)%l)
!        endif
!
!        if(Nr_order>0)then
!            if(Nr==1)then
!                allocate(t_deriv_r_1::Ham%deriv(i_mode)%r)
!            elseif(Nr>1)then
!                allocate(Ham%deriv(i_mode)%r,source=t_deriv_r_N(i_mode))
!            endif
!        else
!            allocate(t_deriv_r_null::Ham%deriv(i_mode)%r)
!        endif
!    enddo
!end subroutine

subroutine set_deriv(Ham)
    class(t_H)      :: Ham(:)
    integer         :: i

    do i=1,size(Ham)
        Call Ham(i)%set_deriv()
    enddo
end subroutine

end module
