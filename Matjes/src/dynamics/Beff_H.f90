module m_Beff_H
use m_derived_types
use m_H_public, only : t_H
implicit none
private
public get_B,get_B_single


contains

subroutine get_B(Ham,lat,B)
!calculates the effective internal magnetic field acting on the magnetization for the dynamics
    class(t_h),intent(in)        ::  Ham(:) !list of Hamiltonians to consider
    class(lattice),intent(in)    ::  lat    !lattice containing current order-parameters 
    real(8),intent(inout)        ::  B(:)

    real(8)      :: tmp(size(B))
    integer      :: Nl,Nr !number of order parameters at left/right side of Hamiltonian 
    integer      :: Nl_m,Nr_m !number magnetic order parameter occurances at left/right side of Hamiltonian 
    integer      :: i

    B=0.0d0
    do i=1,size(Ham)
        Nl=size(Ham(i)%op_l); Nr=size(Ham(i)%op_r)
        Nl_m=count(Ham(i)%op_l==1); Nr_m=count(Ham(i)%op_r==1)
        if(Nl==1.and.Nr==1)then
        !normal rank 2 Hamiltonian
            tmp=0.0d0
            if(Nl_m==1.and.Nr_m==1)then
                Call Ham(i)%mult_r(lat,tmp)
                B=B-2.0d0*tmp
            elseif(Nl_m==1)then
                Call Ham(i)%mult_r(lat,tmp)
                B=B-tmp
            elseif(Nr_m==1)then
                Call Ham(i)%mult_l(lat,tmp)
                B=B-tmp
            endif
        else
        !Hamiltonian with rank higher than 2
            tmp=0.0d0
            !consider left side
            if(Nl_m==1)then
                if(Nl==1)then
                    Call Ham(i)%mult_r(lat,tmp)
                    B=B-tmp
                elseif(Nl_m==1)then
                    Call Ham(i)%mult_r_red(lat,tmp,1) !1 for reduce everything but M
                    B=B-tmp
                else
                    STOP "implement effective magnetic field for Hamiltonian with more than one M occurancy on other side"
                endif
            else
                STOP "implement effective magnetic field for Hamiltonian with more than one M occurancy on a side"
            endif
            !consider right side
            if(Nr_m==1)then
                if(Nr==1)then
                    Call Ham(i)%mult_l(lat,tmp)
                    B=B-tmp
                elseif(Nr_m==1)then
                    Call Ham(i)%mult_l_red(lat,tmp,1) !1 for reduce everything but M
                    B=B-tmp
                else
                    STOP "implement effective magnetic field for Hamiltonian with more than one M occurancy on other side"
                endif
            else
                STOP "implement effective magnetic field for Hamiltonian with more than one M occurancy on a side"
            endif
        endif
    enddo

end subroutine

subroutine get_B_single(Ham,i_site,lat,B)
!calculates the effective internal magnetic field acting on the magnetization for the dynamics
    class(t_h),intent(in)        :: Ham(:) !list of Hamiltonians to consider
    integer,intent(in)           :: i_site
    class(lattice),intent(in)    :: lat    !lattice containing current order-parameters 
    real(8),intent(inout)        :: B(:)

    real(8)      :: tmp(size(B))
    integer      :: Nl,Nr !number of order parameters at left/right side of Hamiltonian 
    integer      :: Nl_m,Nr_m !number magnetic order parameter occurances at left/right side of Hamiltonian 
    integer      :: i

    B=0.0d0
    do i=1,size(Ham)
        Nl=size(Ham(i)%op_l); Nr=size(Ham(i)%op_r)
        Nl_m=count(Ham(i)%op_l==1); Nr_m=count(Ham(i)%op_r==1)
        if(Nl==1.and.Nr==1)then
        !normal rank 2 Hamiltonian
            tmp=0.0d0
            if(Nl_m==1.and.Nr_m==1)then
                Call Ham(i)%mult_r_single(i_site,lat,tmp)
                B=B-2.0d0*tmp
            elseif(Nl_m==1)then
                Call Ham(i)%mult_r_single(i_site,lat,tmp)
                B=B-tmp
            elseif(Nr_m==1)then
                Call Ham(i)%mult_l_single(i_site,lat,tmp)
                B=B-tmp
            endif
        else
        !Hamiltonian with rank higher than 2
            tmp=0.0d0
            !consider left side
            if(Nl_m==1)then
                if(Nl==1)then
                    Call Ham(i)%mult_r_single(i_site,lat,tmp)
                    B=B-tmp
                elseif(Nl_m==1)then
                    Call Ham(i)%mult_r_red_single(i_site,lat,tmp,1) !1 for reduce everything but M
                    B=B-tmp
                else
                    STOP "implement effective magnetic field for Hamiltonian with more than one M occurancy on other side"
                endif
            else
                STOP "implement effective magnetic field for Hamiltonian with more than one M occurancy on a side"
            endif
            !consider right side
            if(Nr_m==1)then
                if(Nr==1)then
                    Call Ham(i)%mult_l_single(i_site,lat,tmp)
                    B=B-tmp
                elseif(Nr_m==1)then
                    Call Ham(i)%mult_l_red_single(i_site,lat,tmp,1) !1 for reduce everything but M
                    B=B-tmp
                else
                    STOP "implement effective magnetic field for Hamiltonian with more than one M occurancy on other side"
                endif
            else
                STOP "implement effective magnetic field for Hamiltonian with more than one M occurancy on a side"
            endif
        endif
    enddo

end subroutine


end module


