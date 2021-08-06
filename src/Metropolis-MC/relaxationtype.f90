module m_relaxtyp
    use m_derived_types, only : lattice
!    use m_H_public
    use m_hamiltonian_collection, only: hamiltonian
    use m_work_ham_single, only:  work_ham_single
    implicit none
    private
    public :: underrelax, overrelax
    
    contains
    ! functions that relaxes the spins with respect of the dE/dM
    ! in one case, the spins are choosen in the direction of -dE/DM so energy diminishes
    ! in the second case, the spins are choosen in the direction of +dE/DM so energy increases
    !
    function underrelax(iomp,lat,work,H)
        !update magnetization according to M=M-Mx(MxB)
        !use identity to avoid cross products : Mx(MxB) = M*(M.B)-B (as M.M=1) ->  M-Mx(MxB)=M*(1-M.B)+B
        use m_vector, only : cross
        implicit none
        ! external variable
        integer, intent(in)                 :: iomp
        type(lattice),intent(in)            :: lat
        type(work_ham_single),intent(inout) :: work
        type(hamiltonian), intent(inout)    :: H
        ! value of the function
        real(8)     :: underrelax(3)
        !internal variable
        real(8)     :: norm_local,dotprod, S_int(3)
        real(8)     :: dumy(3)
    

        Call H%get_eff_field_single(lat,iomp,S_int,work,1,dumy)
        dotprod=DOT_PRODUCT(lat%M%modes_in(:,iomp),S_int)
        underrelax=lat%M%modes_in(:,iomp)*(1.0d0-dotprod)+S_int
        norm_local=norm2(underrelax)
        underrelax=underrelax/norm_local

        !I think underrelax can never become 0 so no check for the norm divison should be fine
        !it would require:
        !    M(1-M.B)+B !=0
        !=>  M+B != M.B*M      !M.B=cos(a)*b with a angle between M and B and b=|B|
        !=>  M+B != b*cos(a)*M   !M+B=x*M can only be possible if B is parallel or antiparallel to M => cos(a)=+- 1; B=b*M
        !=>  (1+b)*M != +- b*M  !which can never be satisfield for |M|=1 for either for cos(a)=1 or cos(a)=-1
        !The only reason to check the norm would be if M is set to zero eg simulating a vacancy
  
        !old implementation (also some norm check was in there)
        !M-Mx(MxB)
        !dumy=cross(lat%M%modes_in(:,iomp),S_int,1,3)
        !S_int=lat%M%modes_in(:,iomp)-cross(lat%M%modes_in(:,iomp),dumy,1,3)
        !norm_local=norm2(S_int)
        !underrelax=S_int/norm_local
    end function underrelax
    
    function overrelax(iomp,lat,work,H)
        !Similar to underrelax, but with M=M+Mx(MxB) -> M=M*(1+M.B)-B
        integer, intent(in)                 :: iomp
        type(lattice),intent(in)            :: lat
        type(work_ham_single),intent(inout) :: work
        type(hamiltonian), intent(inout)    :: H
        real(8)                             :: overrelax(3)
        !internal variable
        real(8)     :: norm_local,dotprod, S_int(3)
        real(8)     :: dumy(3)

        Call H%get_eff_field_single(lat,iomp,S_int,work,1,dumy)
        dotprod=DOT_PRODUCT(lat%M%modes_in(:,iomp),S_int)
        overrelax=lat%M%modes_in(:,iomp)*(1.0d0+dotprod)-S_int
        norm_local=norm2(overrelax)
        overrelax=overrelax/norm_local
    end function overrelax
end module
