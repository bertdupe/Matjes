module m_anisotropy_heisenberg
private
public :: get_anisotropy_H,read_anisotropy_input

contains

subroutine read_anisotropy_input(io_param,fname,io_aniso)
    use m_io_utils
    use m_input_H_types, only: io_H_aniso
    integer,intent(in)              :: io_param
    character(len=*), intent(in)    :: fname
    type(io_H_aniso),intent(out)    :: io_aniso
    !internal
    integer         :: max_entry
    real(8)         :: const_mult

    max_entry=max_ind_variable(io_param,'ani_',fname)
    if(max_entry==0)then
        io_aniso%is_set=.false.
        return
    endif
    io_aniso%is_set=.true.
    allocate(io_aniso%val(3*max_entry),source=0.0d0)
    call get_coeff(io_param,fname,'ani_',io_aniso%val,3)
    const_mult=1.0d0
    call get_parameter(io_param,fname,'c_ani',const_mult)
    io_aniso%val=io_aniso%val*const_mult
end subroutine

subroutine get_anisotropy_H(Ham,io,lat)
    !get anisotropy in t_H Hamiltonian format
    use m_H_public, only: t_H
    use m_derived_types, only: lattice
    use m_setH_util,only: get_coo
    use m_input_H_types, only: io_H_aniso

    class(t_H),intent(inout)    :: Ham
    type(io_H_aniso),intent(in) :: io
    type(lattice),intent(in)    :: lat
    !local 
    integer :: i
    real(8),allocatable :: Htmp(:,:)
    real(8),allocatable :: val_tmp(:)
    integer,allocatable :: ind_tmp(:,:)
    integer,allocatable :: line(:,:)

    if(io%is_set)then
        allocate(line(1,lat%Ncell),source=0) !1, since always onsite
        allocate(Htmp(lat%M%dim_mode,lat%M%dim_mode),source=0.d0)
        do i=1,size(io%val)
            Htmp(i,i)=io%val(i)
        enddo
        !get local Hamiltonian in coo format
        Call get_coo(Htmp,val_tmp,ind_tmp)
        !Anisotropy only has simple onsite terms
        do i=1,lat%Ncell
            line(1,i)=i
        enddo
        Call Ham%init_1(line,val_tmp,ind_tmp,[1,1],lat)
    endif
end subroutine

end module m_anisotropy_heisenberg
