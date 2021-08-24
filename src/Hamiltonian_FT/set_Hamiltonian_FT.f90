module m_set_Hamiltonian_FT
use m_construction_Hk
contains

! routines that prepares the Hamiltonian for the FT
! 1 collect all Hamiltonians of rank 2 and the position  difference
!
subroutine set_Hamiltonians(Ham_res,Ham_comb,keep_res,H_io,lat)
    use m_derived_types
    use m_input_H_types
    use m_get_position,only :get_position_ND_to_1D
    use, intrinsic :: iso_fortran_env, only : output_unit

    use m_symmetry_operators
    use m_anisotropy_heisenberg,only: get_anisotropy_H
    use m_zeeman,only: get_zeeman_H
    use m_exchange_heisenberg_J, only: get_exchange_J
    use m_exchange_heisenberg_D, only: get_exchange_D
    use m_harmonic_phonon,only: get_Forces_F
    class(Hk_inp_t),allocatable,intent(out)  :: Ham_res
    class(Hk_inp_t),allocatable,intent(out)  :: Ham_comb
    logical,intent(in)                  :: keep_res ! keeps the Ham_res terms allocated
    type(io_h),intent(in)               :: H_io
    type(lattice), intent(in)           :: lat

    integer :: i_H,N_ham
    logical :: use_Ham(5)


    use_ham(1)=H_io%J%is_set
    use_ham(2)=H_io%D%is_set
    use_ham(3)=H_io%aniso%is_set
    use_ham(4)=H_io%zeeman%is_set
    use_ham(7)=H_io%F%is_set

    N_ham=count(use_ham)
    allocate(Ham_res(N_Ham),Ham_comb(N_Ham))
!    Call get_Htype_N(Ham_res%H,N_ham)
    i_H=1
    !exchange_J (without DMI)
    if(use_ham(1))then
        Call get_exchange_J(Ham_res(i_H)%H_folded,H_io%J,lat)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif
    !exchange_D (only DMI)
    if(use_ham(2))then
        Call get_exchange_D(Ham_res(i_H)%H_folded,H_io%D,lat)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif
    !anisotropy
    if(use_ham(3))then
        Call get_anisotropy_H(Ham_res(i_H)%H_folded,H_io%aniso,lat)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif
    !zeeman
    if(use_ham(4))then
        Call get_zeeman_H(Ham_res(i_H)%H_folded,H_io%zeeman,lat)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif
    !Harmonic phonon (F)
    if(use_ham(5))then
        Call get_Forces_F(Ham_res(i_H)%H_folded,H_io%F,lat,Ham_res(i_H)%H,Ham_res(i_H)%diffR)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif

    write(output_unit,'(//A,I3,A)') "The following ",N_ham," Hamiltonians in real-space have been set:"
    do i_H=1,N_ham
        write(output_unit,'(3XA)') trim(Ham_res(i_H)%desc)
    enddo
    write(output_unit,'(//)')

    do i_H=1,N_ham
        if(.not. Ham_res(i_h)%is_set()) STOP "not all Hamiltonians are set"
        Call Ham_res(i_h)%optimize()
        Call Ham_res(i_h)%finish_setup()
    enddo

    Call combine_Hamiltonians(keep_res,Ham_res,Ham_comb)
    do i_H=1,size(Ham_comb)
        Call Ham_comb(i_h)%optimize()
        Call Ham_comb(i_h)%finish_setup()
    enddo
end subroutine
end module m_set_Hamiltonian_FT
