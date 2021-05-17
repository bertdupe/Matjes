module m_eff_field
!shall be replaced everywhere by get_eff_field from m_hamiltonian_collection's hamiltonian 
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
    integer,intent(in)           ::  Ham_type   !integer that decides with respect to which mode the Hamiltonians derivative shall be obtained [1,number_different_order_parameters]

    integer     :: iH,NH
    real(8)     :: tmp(size(B))

    B=0.d0
    NH=size(ham)
    do iH=1,NH
        Call Ham(iH)%deriv(Ham_type)%get(Ham(iH),lat,B,tmp)
    enddo
    B=-B    !field is negative derivative
end subroutine

subroutine get_eff_field_single(Ham,i_site,lat,B,Ham_type)
!calculates the effective internal magnetic field acting on the magnetization for the dynamics
    class(t_h),intent(in)        :: Ham(:) !list of Hamiltonians to consider
    integer,intent(in)           :: i_site,Ham_type
    class(lattice),intent(in)    :: lat    !lattice containing current order-parameters 
    real(8),intent(inout)        :: B(:)

    real(8)      :: tmp(size(B))
    integer      :: iH

    B=0.d0
    do iH=1,size(ham)
        Call Ham(iH)%deriv(Ham_type)%get_single(Ham(iH),lat,i_site,B,tmp)
    enddo
    B=-B    !field is negative derivative
end subroutine

end module
