module m_set_Hamiltonian_FT
use m_FT_Ham_coo_rtok_base
use m_H_public
implicit none

contains

! routines that prepares the Hamiltonian for the FT
! 1 collect all Hamiltonians of rank 2 and the position  difference
!
subroutine set_Hamiltonians_FT(FT_Ham,H_io,lat)
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
    use m_exchange_heisenberg_general, only : get_exchange_ExchG

    type(H_inp_real_to_k),allocatable,intent(inout)       :: FT_Ham(:)
    type(io_h),intent(in)                                 :: H_io
    type(lattice), intent(in)                             :: lat

    ! internal
    integer :: i_H,N_ham
    logical :: use_Ham(6)


    use_ham(1)=H_io%J%is_set
    use_ham(2)=H_io%D%is_set
    use_ham(3)=H_io%aniso%is_set
    use_ham(4)=H_io%zeeman%is_set
    use_ham(5)=H_io%F%is_set
    use_ham(6)=H_io%Exchten%is_set

    N_ham=count(use_ham)
    if (.not.allocated(FT_Ham)) allocate(FT_Ham(N_ham))

    i_H=1
    !exchange_J (without DMI)
    if(use_ham(1))then
        call get_Htype(FT_Ham(i_H)%H_folded)
        Call get_exchange_J(FT_Ham(i_H)%H_folded,H_io%J,lat,FT_Ham(i_H)%H,FT_Ham(i_H)%diffR)
        if(FT_Ham(i_H)%H_folded%is_set()) i_H=i_H+1
    endif
    !exchange_D (only DMI)
    if(use_ham(2))then
        call get_Htype(FT_Ham(i_H)%H_folded)
        Call get_exchange_D(FT_Ham(i_H)%H_folded,H_io%D,lat,FT_Ham(i_H)%H,FT_Ham(i_H)%diffR)
        if(FT_Ham(i_H)%H_folded%is_set()) i_H=i_H+1
    endif
    !anisotropy
    if(use_ham(3))then
        call get_Htype(FT_Ham(i_H)%H_folded)
        Call get_anisotropy_H(FT_Ham(i_H)%H_folded,H_io%aniso,lat,FT_Ham(i_H)%H,FT_Ham(i_H)%diffR)
        if(FT_Ham(i_H)%H_folded%is_set()) i_H=i_H+1
    endif
    !zeeman
    if(use_ham(4))then
        call get_Htype(FT_Ham(i_H)%H_folded)
        Call get_zeeman_H(FT_Ham(i_H)%H_folded,H_io%zeeman,lat,FT_Ham(i_H)%H,FT_Ham(i_H)%diffR)
        if(FT_Ham(i_H)%H_folded%is_set()) i_H=i_H+1
    endif
    !Harmonic phonon (F)
    if(use_ham(5))then
        call get_Htype(FT_Ham(i_H)%H_folded)
        Call get_Forces_F(FT_Ham(i_H)%H_folded,H_io%F,lat,FT_Ham(i_H)%H,FT_Ham(i_H)%diffR)
        if(FT_Ham(i_H)%H_folded%is_set()) i_H=i_H+1
    endif
    !General exchange tensor
    if(use_ham(6))then
        call get_Htype(FT_Ham(i_H)%H_folded)
        Call get_exchange_ExchG(FT_Ham(i_H)%H_folded,H_io%Exchten,lat,FT_Ham(i_H)%H,FT_Ham(i_H)%diffR)
        if(FT_Ham(i_H)%H_folded%is_set()) i_H=i_H+1
    endif

    write(output_unit,'(//A,I3,A)') "The following ",N_ham," Hamiltonians FT have been set:"
    do i_H=1,N_ham
        write(output_unit,'(3XA)') trim(FT_Ham(i_H)%H_folded%desc)
    enddo
    write(output_unit,'(//)')

!    Call combine_Hamiltonians(keep_res,Ham_res,Ham_comb)
!    do i_H=1,size(Ham_comb)
!        Call Ham_comb(i_h)%optimize()
!        Call Ham_comb(i_h)%finish_setup()
!    enddo

end subroutine

end module m_set_Hamiltonian_FT
