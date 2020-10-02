module m_Beff_H
use m_derived_types
use m_H_type
implicit none


contains

subroutine get_B(Ham,lat,B)
!calculates the effective internal magnetic field acting on the magnetization for the dynamics
    class(t_h),intent(in)        ::  Ham(:) !list of Hamiltonians to consider
    class(lattice),intent(in)    ::  lat    !lattice containing current order-parameters 
    real(8),intent(inout)        ::  B(:)

    real(8)      :: tmp(size(B))
    integer      :: NH
    integer      :: Nl,Nr
    integer      :: i

    NH=size(Ham)
    B=0.0d0
    do i=1,NH
        if(size(Ham(i)%op_l)>1.or.size(Ham(i)%op_r)>1) STOP "implement more than rank2 operator for get_B in m_Beff_H"

        Nl=count(Ham(i)%op_l==1);Nr=count(Ham(i)%op_r==1)
        tmp=0.0d0
        if(Nl==1.and.Nr==1)then
            Call Ham(i)%mult_r(lat,tmp)
            tmp=-2.0d0*tmp
        elseif(Nl==1)then
            Call Ham(i)%mult_r(lat,tmp)
            tmp=-tmp
        elseif(Nr==1)then
            Call Ham(i)%mult_l(lat,tmp)
            tmp=-tmp
        else
            continue
        endif
        B=B+tmp
    enddo

end subroutine

end module


