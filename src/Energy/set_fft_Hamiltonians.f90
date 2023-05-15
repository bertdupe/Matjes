module m_set_fft_Hamiltonians
use m_fft_H_public
implicit none
private
public :: set_fft_Hamiltonians
contains

subroutine set_fft_Hamiltonians(Ham_res,Ham_comb,keep_res,H_io,lat)
    use m_input_H_types 
    use m_derived_types
    use m_H_type, only: len_desc
    use, intrinsic :: iso_fortran_env, only : output_unit

    use m_exchange_heisenberg_J, only: get_exchange_J_fft
    use m_exchange_heisenberg_D, only: get_exchange_D_fft
    use m_anisotropy_heisenberg, only: get_anisotropy_fft
    use m_dipolar_magnetic, only: get_dipolar_fft
    use m_dipolar_phonon, only: get_dipolar_ph_fft
    use m_exchange_heisenberg_general, only : get_exchange_ExchG_fft

    class(fft_H),allocatable,intent(out)    :: Ham_res(:)
    class(fft_H),allocatable,intent(out)    :: Ham_comb(:)
    logical,intent(in)                      :: keep_res ! keeps the Ham_res terms allocated
    type(io_h),intent(in)                   :: H_io
    type(lattice), intent(in)               :: lat
    character(len=len_desc)     :: desc=""
    integer                     :: i_H,N_ham
    logical                     :: use_Ham(7)

    use_ham(1)=H_io%J%is_set.and.H_io%J%fft
    use_ham(2)=H_io%D%is_set.and.H_io%D%fft
    use_ham(3)=H_io%aniso%is_set.and.H_io%aniso%fft
    use_ham(4)=H_io%dip%is_set.and.H_io%dip%fft
    use_ham(5)=H_io%Exchten%is_set.and.H_io%Exchten%fft
    use_ham(6)=H_io%SC%is_set.and.H_io%SC%fft
    use_ham(7)=H_io%dip_ph%is_set.and.H_io%dip_ph%fft

    N_ham=count(use_ham)
    if(N_ham<1) return  !nothing to do here
    Call get_fft_H_N(Ham_res,N_ham)
    i_H=1 
    !exchange_J (without DMI)
    if(use_ham(1))then
        Call get_exchange_J_fft(Ham_res(i_H),H_io%J,lat)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif
    !exchange_D (only DMI)
    if(use_ham(2))then
        Call get_exchange_D_fft(Ham_res(i_H),H_io%D,lat)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif
    !anisotropy term interaction
    if(use_ham(3))then
        Call get_anisotropy_fft(Ham_res(i_H),H_io%aniso,lat)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif
    !magnetic dipolar interaction
    if(use_ham(4))then
        Call get_dipolar_fft(Ham_res(i_H),H_io%dip,lat)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif
    !general exchange tensor
    if(use_ham(5))then
        Call get_exchange_ExchG_fft(Ham_res(i_H),H_io%Exchten,lat)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif
    
    !electric dipolar interaction
    if(use_ham(7))then
        Call get_dipolar_ph_fft(Ham_res(i_H),H_io%dip_ph,lat)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif
    
    !spin current DMI
!    if(use_ham(6))then
!        Call get_coupling_SC_fft(Ham_res(i_H),H_io%SC,lat)
!        if(Ham_res(i_H)%is_set()) i_H=i_H+1
!    endif

    write(output_unit,'(//A,I3,A)') "The following ",N_ham," Hamiltonians in fourier-space have been set:"
    do i_H=1,N_ham
        Call Ham_res(i_H)%get_desc(desc)
        write(output_unit,'(3XA)') trim(desc)
    enddo
    write(output_unit,'(//)')

    Call combine_Hamiltonians(keep_res,Ham_res,Ham_comb)
end subroutine

subroutine combine_Hamiltonians(keep_res,Ham_res,Ham_comb)
    !combines the fft Hamiltonians assuming that they all act on the same space
    class(fft_H),allocatable,intent(inout)    :: Ham_res(:)
    class(fft_H),allocatable,intent(inout)    :: Ham_comb(:)
    logical,intent(in)                      :: keep_res ! keeps the Ham_res terms allocated

    integer     ::  i

    Call get_fft_H_N(Ham_comb,1)
    do i=1,size(ham_res)
        Call Ham_comb(1)%add(Ham_res(i))
    enddo
    if(.not.keep_res)then
        do i=1,size(ham_res)
            Call Ham_res(i)%destroy()
        enddo
        deallocate(Ham_res)
    endif
end subroutine
end module
