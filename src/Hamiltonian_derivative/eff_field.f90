module m_eff_field
use m_derived_types
use m_H_public, only : t_H
implicit none
private
public :: get_eff_field,get_eff_field_single


contains

subroutine get_eff_field(Ham,lat,B,Ham_type)
!calculates the effective internal magnetic field acting on the magnetization for the dynamics
    class(t_h),intent(in)        ::  Ham(:) !list of Hamiltonians to consider
    class(lattice),intent(in)    ::  lat    !lattice containing current order-parameters 
    real(8),intent(inout)        ::  B(:)
    integer,intent(in)           ::  Ham_type   !integer that decides which Hamiltonian should be derived

    real(8)      :: tmp(size(B))
    integer      :: Nl,Nr !number of order parameters at left/right side of Hamiltonian 
    integer      :: Nl_order,Nr_order !numbern of order parameters occurances at left/right side of Hamiltonian
    integer      :: i

    B=0.0d0
    do i=1,size(Ham)
        Nl=size(Ham(i)%op_l); Nr=size(Ham(i)%op_r)
        Nl_order=count(Ham(i)%op_l==Ham_type); Nr_order=count(Ham(i)%op_r==Ham_type)
        if(Nl==1.and.Nr==1)then
        !normal rank 2 Hamiltonian
            tmp=0.0d0
            if(Nl_order==1.and.Nr_order==1)then
                Call Ham(i)%mult_r(lat,tmp)
                B=B-2.0d0*tmp
            elseif(Nl_order==1)then
                Call Ham(i)%mult_r(lat,tmp)
                B=B-tmp
            elseif(Nr_order==1)then
                Call Ham(i)%mult_l(lat,tmp)
                B=B-tmp
            endif
        else
        !Hamiltonian with rank higher than 2
            tmp=0.0d0
            !consider left side
            if(Nl_order==1)then
                if(Nl==1)then
                    Call Ham(i)%mult_r(lat,tmp)
                    B=B-tmp
                elseif(Nl_order==1)then
                    Call Ham(i)%mult_r_red(lat,tmp,Ham_type) !1 for reduce everything but ham_type 
                    B=B-tmp
                else
                    STOP "implement effective field for Hamiltonian with more than one occurance of Ham_type on other side"
                endif
            else
                STOP "implement effective field for Hamiltonian with more than one occurance of Ham_type on a side"
            endif
            !consider right side
            if(Nr_order==1)then
                if(Nr==1)then
                    Call Ham(i)%mult_l(lat,tmp)
                    B=B-tmp
                elseif(Nr_order==1)then
                    Call Ham(i)%mult_l_red(lat,tmp,Ham_type) !1 for reduce everything but Ham_type 
                    B=B-tmp
                else
                    STOP "implement effective field for Hamiltonian with more than one occurance of Ham_type on other side"
                endif
            else
                STOP "implement effective field for Hamiltonian with more than one occurance of Ham_type on a side"
            endif
        endif
    enddo

end subroutine

subroutine get_eff_field_single(Ham,i_site,lat,B,Ham_type)
!calculates the effective internal magnetic field acting on the magnetization for the dynamics
    class(t_h),intent(in)        :: Ham(:) !list of Hamiltonians to consider
    integer,intent(in)           :: i_site,Ham_type
    class(lattice),intent(in)    :: lat    !lattice containing current order-parameters 
    real(8),intent(inout)        :: B(:)

    real(8)      :: tmp(size(B))
    integer      :: Nl,Nr !number of order parameters at left/right side of Hamiltonian 
    integer      :: Nl_order,Nr_order !number of order parameters occurances at left/right side of Hamiltonian
    integer      :: i

    B=0.0d0
    do i=1,size(Ham)
        Nl=size(Ham(i)%op_l); Nr=size(Ham(i)%op_r)
        Nl_order=count(Ham(i)%op_l==Ham_type); Nr_order=count(Ham(i)%op_r==Ham_type)
        if(Nl==1.and.Nr==1)then
        !normal rank 2 Hamiltonian
            tmp=0.0d0
            if(Nl_order==1.and.Nr_order==1)then
                Call Ham(i)%mult_r_single(i_site,lat,tmp)
                B=B-2.0d0*tmp
            elseif(Nl_order==1)then
                Call Ham(i)%mult_r_single(i_site,lat,tmp)
                B=B-tmp
            elseif(Nr_order==1)then
                Call Ham(i)%mult_l_single(i_site,lat,tmp)
                B=B-tmp
            endif
        else
        !Hamiltonian with rank higher than 2
            tmp=0.0d0
            !consider left side
            if(Nl_order==1)then
                if(Nl==1)then
                    Call Ham(i)%mult_r_single(i_site,lat,tmp)
                    B=B-tmp
                elseif(Nl_order==1)then
                    Call Ham(i)%mult_r_red_single(i_site,lat,tmp,Ham_type) !reduce everything but ham_type 
                    B=B-tmp
                else
                    STOP "implement effective field for Hamiltonian with more than one occurance of Ham_type on other side"
                endif
            else
                STOP "implement effective field for Hamiltonian with more than one occurance of Ham_type on a side"
            endif
            !consider right side
            if(Nr_order==1)then
                if(Nr==1)then
                    Call Ham(i)%mult_l_single(i_site,lat,tmp)
                    B=B-tmp
                elseif(Nr_order==1)then
                    Call Ham(i)%mult_l_red_single(i_site,lat,tmp,Ham_type) !reduce everything but ham_type 
                    B=B-tmp
                else
                    STOP "implement effective field for Hamiltonian with more than one occurance of Ham_type on other side"
                endif
            else
                STOP "implement effective field for Hamiltonian with more than one occurance of Ham_type on a side"
            endif
        endif
    enddo

end subroutine


end module


