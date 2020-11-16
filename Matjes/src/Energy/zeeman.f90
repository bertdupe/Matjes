module m_zeeman
implicit none
private
public :: get_zeeman_H,read_zeeman_input

contains

subroutine read_zeeman_input(io_param,fname,io)
    use m_io_utils
    use m_input_H_types, only: io_H_zeeman
    integer,intent(in)              :: io_param
    character(len=*), intent(in)    :: fname
    type(io_H_zeeman),intent(out)   :: io
    !internal
    real(8)         :: h_ext(3)
    logical         :: enable_zeeman

    call get_parameter(io_param,fname,'c_zeeman',io%c_zeeman)
    h_ext=0.d0
    call get_parameter(io_param,fname,'H_ext',3,h_ext)
    enable_zeeman=.false.
    call get_parameter(io_param,fname,'enable_zeeman',enable_zeeman)
    io%is_set=enable_zeeman.or.norm2(h_ext).ge.1.0d-8
end subroutine

subroutine get_zeeman_H(Ham,io,lat)
    !get zeeman in t_H Hamiltonian format
    use m_H_public
    use m_derived_types
    use m_setH_util,only: get_coo
    use m_input_H_types, only: io_H_zeeman
    use m_constants, only : mu_B,mu_0

    class(t_H),intent(inout)    :: Ham
    type(io_H_zeeman),intent(in):: io
    type(lattice),intent(in)    :: lat
    !local parameters
    real(8),allocatable     :: magmom(:)
    integer                 :: i
    !describe local Hamiltonian
    real(8),allocatable     :: Htmp(:,:)
    !input for setting t_H
    real(8),allocatable     :: val_tmp(:)
    integer,allocatable     :: ind_tmp(:,:)
    integer,allocatable     :: line(:,:)

    if(io%is_set)then
        allocate(line(1,lat%Ncell),source=0) !1, since always onsite
        allocate(Htmp(3,lat%M%dim_mode),source=0.d0) !assume shape of B-field has to be 3
        Call lat%cell%get_mag_magmom(magmom)
        do i=1,size(magmom)
            Htmp(1,(i-1)*3+1)=magmom(i)
            Htmp(2,(i-1)*3+2)=magmom(i)
            Htmp(3,(i-1)*3+3)=magmom(i)
        enddo
        Htmp=(mu_0*mu_B*io%c_zeeman) * Htmp

        !get local Hamiltonian in coo format
        Call get_coo(Htmp,val_tmp,ind_tmp)
        !simple on-site term
        do i=1,lat%Ncell
            line(1,i)=i
        enddo
        Call Ham%init_1(line,val_tmp,ind_tmp,[3,1],lat)
        Ham%desc="Zeeman energy"
    endif
end subroutine

end module m_zeeman
