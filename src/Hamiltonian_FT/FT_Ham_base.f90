module m_FT_Ham_base
use m_H_combined, only : t_h
implicit none

private
public :: H_inp_real_to_k

type H_inp_real_to_k
! folded Hamiltonian that contains the information about the bounds (useless in the case of the FT)
    class(t_h),allocatable       :: H_folded
! Hamiltonian resolved in shell and bound
    class(t_h),allocatable       :: H(:)
    real(8),allocatable         :: diffR(:,:)
contains
!    procedure :: combine => Hk_inp_combine  !combine 2 Hk_inp_t by copying
!    procedure :: destroy => Hk_inp_destroy
end type

end module
