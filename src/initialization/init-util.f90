module m_init_util
implicit none
public

contains

subroutine get_pos(lat,dim_mode,ordname,pos,dim_state)
    use m_derived_types,only: lattice 
    type(lattice), intent(in)       :: lat      !entire lattice containing geometric information
    integer,intent(in)              :: dim_mode
    character(*),intent(in)         :: ordname  !name of the order parameter
    real(8),allocatable,intent(out) :: pos(:)
    integer,intent(out)             :: dim_state
    integer     :: nmag
    
    nmag=lat%cell%num_mag()
    if(dim_mode==1.and.nmag==1)then
        !choose position of magnetic atoms for initialization
        Call lat%get_pos_center(pos)
        dim_state=1
    elseif(dim_mode==1.and.nmag>1)then
        !use center of cell as position
        Call lat%get_pos_mag(pos)
        dim_state=1
    elseif(dim_mode==nmag*3)then
        !choose position of magnetic atoms for initialization
        Call lat%get_pos_mag(pos)
        dim_state=3
    elseif(dim_mode==3)then
        !use center of cell as position
        Call lat%get_pos_center(pos)
        dim_state=3
        write(*,*) "CANNOT USE get_pos FOR",ordname,"as its dim_modes is most unexpected" !add something with positions of all atoms, and multiplied by 3? (risky mixup in natom=3*nmag
        return
    endif
end subroutine

subroutine get_pos_vec(lat,dim_mode,ordname,pos)
    use m_derived_types,only: lattice 
    type(lattice), intent(in)       :: lat      !entire lattice containing geometric information
    integer,intent(in)              :: dim_mode
    character(*),intent(in)         :: ordname  !name of the order parameter
    real(8),allocatable,intent(out) :: pos(:)
    integer     :: nmag
    
    nmag=lat%cell%num_mag()
    if(dim_mode==nmag*3)then
        !choose position of magnetic atoms for initialization
        Call lat%get_pos_mag(pos)
    elseif(dim_mode==3)then
        !use center of cell as position
        Call lat%get_pos_center(pos)
    else
        write(*,*) "CANNOT USE GET_POS_VEC FOR",ordname,"as its dim_mode is neither 3 nor 3*nmag"
        return
    endif
end subroutine

end module

