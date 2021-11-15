module m_FT_Ham_coo_rtok_base
use m_H_combined, only : t_h
use m_H_type
implicit none

private
public H_inp_real_to_k

type  H_inp_real_to_k
! folded Hamiltonian that contains the information about the bounds (useless in the case of the FT)
    class(t_h),allocatable      :: H_folded
! Hamiltonian resolved in shell and bound
    real(8),allocatable         :: H(:,:,:)
    real(8),allocatable         :: diffR(:,:)
contains
!    procedure :: combine => Hk_inp_combine  !combine 2 Hk_inp_t by copying
    procedure   :: destroy => H_inp_real_to_k_destroy
end type

contains

subroutine H_inp_real_to_k_destroy(this)
    class(H_inp_real_to_k),intent(inout)   ::  this
    integer ::  i

    Call this%H_folded%destroy()
    deallocate(this%H)
    deallocate(this%diffR)
end subroutine

end module
