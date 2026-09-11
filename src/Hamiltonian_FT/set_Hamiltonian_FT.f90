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
    use m_spincurrent, only :get_coupling_SC
    use m_ASR_phonon, only : get_ASR_Ph

    type(H_inp_real_to_k),allocatable,intent(inout)       :: FT_Ham(:)
    class(t_H),allocatable                                :: Ham_comb(:)
    type(io_h),intent(in)                                 :: H_io
    type(lattice), intent(in)                             :: lat

    ! internal
    integer :: i_H,N_ham
    logical :: use_Ham(8)


    use_ham(1)=H_io%J%is_set
    use_ham(2)=H_io%D%is_set
    use_ham(3)=H_io%aniso%is_set
    use_ham(4)=H_io%zeeman%is_set
    use_ham(5)=H_io%F%is_set
    use_ham(6)=H_io%ASR_ph%is_set
    use_ham(7)=H_io%Exchten%is_set
    use_ham(8)=H_io%SC%is_set

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
    !ASR of Harmonic phonon (F)
    if(use_ham(6))then
        call get_Htype(FT_Ham(i_H)%H_folded)
        Call get_ASR_Ph(FT_Ham(i_H)%H_folded,H_io%ASR_Ph,lat,FT_Ham(i_H)%H,FT_Ham(i_H)%diffR)
        if(FT_Ham(i_H)%H_folded%is_set()) i_H=i_H+1
    endif
    !General exchange tensor
    if(use_ham(7))then
        call get_Htype(FT_Ham(i_H)%H_folded)
        Call get_exchange_ExchG(FT_Ham(i_H)%H_folded,H_io%Exchten,lat,FT_Ham(i_H)%H,FT_Ham(i_H)%diffR)
        if(FT_Ham(i_H)%H_folded%is_set()) i_H=i_H+1
    endif
    !spin current model
    if(use_ham(8))then
        call get_Htype(FT_Ham(i_H)%H_folded)
        Call get_coupling_SC(FT_Ham(i_H)%H_folded,H_io%SC,lat,FT_Ham(i_H)%H,FT_Ham(i_H)%diffR)
        if(FT_Ham(i_H)%H_folded%is_set()) i_H=i_H+1
    endif

    write(output_unit,'(//A,I3,A)') "The following ",N_ham," Hamiltonians FT have been set:"
    do i_H=1,N_ham
        write(output_unit,'(3XA)') trim(FT_Ham(i_H)%H_folded%desc)
    enddo
    write(output_unit,'(//)')

    do i_H=1,N_ham
        if(.not. FT_Ham(i_h)%H_folded%is_set()) STOP "not all Hamiltonians are set"
        Call FT_Ham(i_H)%H_folded%optimize()
        Call FT_Ham(i_H)%H_folded%finish_setup()
    enddo

    Call combine_Hamiltonians(FT_Ham,Ham_comb)
    do i_H=1,size(Ham_comb)
        Call Ham_comb(i_h)%optimize()
        Call Ham_comb(i_h)%finish_setup()
    enddo

end subroutine

subroutine combine_Hamiltonians(H_in,H_comb)
    !subroutine which combines Hamiltonians from the input array, if they act on the same space
    !THIS NEEDS TO BE UPDATED IF THE SPACE IS NOT UNIQUE ANYMORE (eg. terms in one dimension that are not onsite)
    class(H_inp_real_to_k),intent(inout)    :: H_in(:)
    class(t_H),allocatable,intent(inout)    :: H_comb(:)
    !internal
    integer         :: N_in,N_out
    integer         :: max_dim
    integer,allocatable ::  space(:,:,:)    !save all space indicators
    integer         ::  i_unique(size(H_in))
    integer         ::  i,j

    !find out which operators are in which space
    N_in=size(H_in)
    max_dim=0
    do i=1,N_in
        max_dim=max(max_dim,size(H_in(i)%H_folded%op_l))
        max_dim=max(max_dim,size(H_in(i)%H_folded%op_r))
    enddo

    allocate(space(max_dim,2,N_in),source=0)
    do i=1,N_in
        do j=1,size(H_in(i)%H_folded%op_l)
            space(j,1,i)=H_in(i)%H_folded%op_l(j)
        enddo
        do j=1,size(H_in(i)%H_folded%op_r)
            space(j,2,i)=H_in(i)%H_folded%op_r(j)
        enddo
    enddo

    N_out=0
    outer:do j=1,N_in
        do i=1,j-1
            if(all(space(:,:,j)==space(:,:,i)))then
                if(H_in(j)%H_folded%same_space(H_in(i)%H_folded))then
                    i_unique(j)=i_unique(i)
                    cycle outer
                endif
            endif
        enddo
        N_out=N_out+1
        i_unique(j)=N_out
    enddo outer

    Call get_Htype_N(H_comb,N_out)
    do i=1,N_in
        Call H_comb(i_unique(i))%add(H_in(i)%H_folded)
    enddo

end subroutine

end module m_set_Hamiltonian_FT
