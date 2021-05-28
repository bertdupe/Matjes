module m_deriv_public
use m_deriv_base, only: t_deriv
use m_H_type, only: t_H_base
use m_derived_types,only : lattice
use m_work_ham_single, only:  work_ham_single, work_mode
use m_deriv_rank2
use m_deriv_rankN
implicit none
private
public t_deriv, set_deriv
contains

subroutine set_deriv(deriv,order,op_l,op_r)
    !sets all the derivative pointers and orders
    type(t_deriv),intent(inout)         :: deriv
    integer,intent(in)                  :: order            !order parameter index
    integer,intent(in)                  :: op_l(:),op_r(:)  !indices of left and right order parameters

    integer         :: Nl,Nr !number of order parameters at left/right side of Hamiltonian 
    integer         :: Nl_order,Nr_order !numbern of order parameters occurances at left/right side of Hamiltonian

    deriv%order=order

    Nl=size(op_l); Nr=size(op_r)
    Nl_order=count(op_l==order); Nr_order=count(op_r==order)

    !set global derivatives
    deriv%l=>get_null
    deriv%r=>get_null
    if(Nl==1.and.Nr==1)then
        !rank 2 Hamiltonian
        if(Nl_order>0) deriv%l=>get_l1
        if(Nr_order>0) deriv%r=>get_r1
        if(Nl_order>0.and.Nr_order>0)then
            !operator has to be symmetric (use only one side)
            deriv%l=>get_l1_sym
            deriv%r=>get_null
        endif
    elseif(Nl>1.and.Nr>1)then
        if(Nl_order>0) deriv%l=>get_lN
        if(Nr_order>0) deriv%r=>get_rN
    elseif(Nl>1.and.Nr==1)then
        if(Nl_order>0) deriv%l=>get_lN
        if(Nr_order>0) deriv%r=>get_r1
    elseif(Nl==1.and.Nr>1)then
        if(Nl_order>0) deriv%l=>get_l1
        if(Nr_order>0) deriv%r=>get_rN
    else
        ERROR STOP "ERROR IN SET_DERIV_NEW, this should never be reached"
    endif

    !set single derivatives
    deriv%l_single=>get_null_single
    deriv%r_single=>get_null_single
    if(Nl==1.and.Nr==1)then
        !rank 2 Hamiltonian
        if(Nl_order>0) deriv%l_single=>get_l1_single
        if(Nr_order>0) deriv%r_single=>get_r1_single
        if(Nl_order>0.and.Nr_order>0)then
            !operator has to be symmetric (use only one side)
            deriv%l_single=>get_l1_sym_single
            deriv%r_single=>get_null_single
        endif
    elseif(Nl>1.and.Nr>1)then
        if(Nl_order>0) deriv%l_single=>get_lN_single
        if(Nr_order>0) deriv%r_single=>get_rN_single
    elseif(Nl>1.and.Nr==1)then
        if(Nl_order>0) deriv%l_single=>get_lN_single 
        if(Nr_order>0) deriv%r_single=>get_r1_single
    elseif(Nl==1.and.Nr>1)then
        if(Nl_order>0) deriv%l_single=>get_l1_single
        if(Nr_order>0) deriv%r_single=>get_rN_single
    else
        ERROR STOP "ERROR IN SET_DERIV_NEW, this should never be reached"
    endif
end subroutine

subroutine get_null(this,H,lat,vec,work)
    class(t_deriv),intent(in)       :: this
    class(t_H_base),intent(in)      :: H
    type(lattice),intent(in)        :: lat
    real(8),intent(inout)           :: vec(:)
    type(work_mode),intent(inout)   :: work
    
    continue
end subroutine

subroutine get_null_single(this,H,lat,site,work,vec)
    class(t_deriv),intent(in)           :: this
    class(t_H_base),intent(in)          :: H
    type(lattice),intent(in)            :: lat
    integer,intent(in)                  :: site
    type(work_ham_single),intent(inout) :: work    !data type containing the temporary data for this calculation to prevent constant allocations/deallocations
    real(8),intent(inout)               :: vec(:)
    
    vec=0.0d0 !OR CONTINUE IF SINGLE NOW WORKS ADDITIVE
    continue    !nothing to do here
end subroutine

end module
