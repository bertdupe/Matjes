module m_relaxtyp
    use m_derived_types, only : lattice
    use m_eff_field, only : get_eff_field_single
    use m_H_public
    implicit none
    private
    public :: underrelax, overrelax
    
    contains
    ! functions that relaxes the spins with respect of the dE/dM
    ! in one case, the spins are choosen in the direction of -dE/DM so energy diminishes
    ! in the second case, the spins are choosen in the direction of +dE/DM so energy increases
    !
    function underrelax(iomp,lat,Hams)
        use m_vector, only : cross
        implicit none
        ! external variable
        integer, intent(in) :: iomp
        type(lattice),intent(in)    :: lat
        class(t_H), intent(in) :: Hams(:)
        ! value of the function
        real(kind=8), dimension(3) :: underrelax
        !internal variable
        real(kind=8), dimension(3) ::S_int
        real(kind=8) :: norm_local,dumy(3)
    
        if(lat%M%dim_mode>3) ERROR STOP "THIS WILL FAIL FOR MORE THAN 1 MAGNETIC MOMENT IN UNITCELL"
        Call get_eff_field_single(Hams,iomp,lat,S_int,1) !1 for magnetism (order_parameter_abbrev)
        norm_local=norm2(S_int)
        
        if (norm_local.gt.1.0d-8) then
           dumy=cross(lat%M%modes_in(:,iomp),S_int,1,3)
           S_int=lat%M%modes_v(:,iomp)-cross(lat%M%modes_in(:,iomp),dumy,1,3)
           norm_local=norm2(S_int)
           underrelax=S_int/norm_local
        else
           underrelax=lat%M%modes_in(:,iomp)
        endif
        
    end function underrelax
    
    function overrelax(iomp,lat,Hams)
        use m_vector, only : cross,norm
        implicit none
        ! external variable
        integer, intent(in) :: iomp
        type(lattice),intent(in)    :: lat
        class(t_H), intent(in) :: Hams(:)
        real(8), dimension(3) :: overrelax
    
        overrelax=-underrelax(iomp,lat,Hams)
        
    end function overrelax
end module
