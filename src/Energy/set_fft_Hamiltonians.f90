module m_set_fft_Hamiltonians
use m_H_public
implicit none
private
public :: set_fft_Hamiltonians
contains

subroutine set_fft_Hamiltonians(Ham_res,Ham_comb,keep_res,H_io,lat)
    use m_input_H_types 
    use m_derived_types

    use m_dipolar_fft
    class(fft_H),allocatable,intent(out)    :: Ham_res(:)
    class(fft_H),allocatable,intent(out)    :: Ham_comb(:)
    logical,intent(in)                      :: keep_res ! keeps the Ham_res terms allocated
    type(io_h),intent(in)                   :: H_io
    type(lattice), intent(in)               :: lat

    integer :: i_H,N_ham
    logical :: use_Ham(1)

    use_ham(1)=H_io%dip%is_set.and.H_io%dip%fft

    N_ham=count(use_ham)
    allocate(fft_H::Ham_res(N_ham))
    i_H=1 
    !dipolar interaction
    if(use_ham(1))then
        Call get_dipolar_fft(Ham_res(i_H),H_io%dip,lat)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif

    Call combine_Hamiltonians(keep_res,Ham_res,Ham_comb)
end subroutine

subroutine combine_Hamiltonians(keep_res,Ham_res,Ham_comb)
    class(fft_H),allocatable,intent(inout)    :: Ham_res(:)
    class(fft_H),allocatable,intent(inout)    :: Ham_comb(:)
    logical,intent(in)                      :: keep_res ! keeps the Ham_res terms allocated

    if(size(Ham_res)/=1) ERROR STOP "ACTUALLY IMPLEMENT THIS"
    if(keep_res)then
        allocate(fft_H::Ham_comb(size(ham_res)))
        Call Ham_res(1)%copy(Ham_comb(1))
    else
        Call move_alloc(Ham_res,Ham_comb)
    endif
end subroutine
end module
