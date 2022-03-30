module m_stark
implicit none
character(len=*),parameter  :: ham_desc="Stark energy"
private
public :: get_stark_E,read_stark_input

contains

subroutine read_stark_input(io_param,fname,io)
    use m_io_utils
    use m_input_H_types, only: io_H_stark
    integer,intent(in)              :: io_param
    character(len=*), intent(in)    :: fname
    type(io_H_stark), intent(out)   :: io
    !internal
    real(8)         :: E_ext(3)
    logical         :: enable_stark

    call get_parameter(io_param,fname,'c_stark',io%c_stark)
    e_ext=0.d0
    call get_parameter(io_param,fname,'E_ext',3,e_ext)
    enable_stark=.false.
    call get_parameter(io_param,fname,'enable_stark',enable_stark)
    io%is_set=enable_stark.or.norm2(e_ext).ge.1.0d-8
end subroutine

subroutine get_stark_E(Ham,io,lat)
    !get stark in t_H Hamiltonian format
    use m_H_public
    use m_derived_types
    use m_setH_util,only : get_coo
    use m_constants, only : J_to_eV
    use m_input_H_types, only : io_H_stark
    use m_constants, only : mu_B,epsilon_0
    use m_mode_public

    class(t_H),intent(inout)    :: Ham
    type(io_H_stark),intent(in) :: io
    type(lattice),intent(in)    :: lat
    !local parameters
    real(8),allocatable     :: Zborn(:)
    integer                 :: i
    !describe local Hamiltonian
    real(8),allocatable     :: Htmp(:,:)
    !input for setting t_H
    real(8),allocatable     :: val_tmp(:)
    integer,allocatable     :: ind_tmp(:,:)
    integer,allocatable     :: connect(:,:)

    if(io%is_set)then
        allocate(Htmp(lat%E%dim_mode,lat%u%dim_mode),source=0.d0) !assume shape of E-field has to be 3
        Call lat%cell%get_Z_phonon(Zborn)
        do i=1,size(Zborn)
            Htmp((i-1)*3+1,(i-1)*3+1)=Zborn(i)
            Htmp((i-1)*3+2,(i-1)*3+2)=Zborn(i)
            Htmp((i-1)*3+3,(i-1)*3+3)=Zborn(i)
        enddo
        Htmp=(io%c_stark) * Htmp * J_to_eV

        !get local Hamiltonian in coo format
        Call get_coo(Htmp,val_tmp,ind_tmp)

        allocate(connect(2,lat%Ncell))
        do i=1,lat%Ncell
            connect(:,i)=i
        enddo
        Call Ham%init_connect(connect,val_tmp,ind_tmp,"EU",lat,1)
        Ham%desc=ham_desc
        !set modes
        Call mode_set_rank1(Ham%mode_l,lat,"E")
        Call mode_set_rank1(Ham%mode_r,lat,"U")
    endif
end subroutine

end module m_stark
