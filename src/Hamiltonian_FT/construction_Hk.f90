module m_construction_Hk
use m_H_public, only : t_h

type Hk_inp_t
! folded Hamiltonian that contains the information about the bounds (useless in the case of the FT)
    class(t_h),allocatable       :: H_folded
! Hamiltonian resolved in shell and bound
    class(t_h),allocatable       :: H(:,:)
    real(8),allocatable         :: diffR(:,:)
contains
!    procedure :: combine => Hk_inp_combine  !combine 2 Hk_inp_t by copying
!    procedure :: destroy => Hk_inp_destroy
end type

contains

!subroutine get_Hk(Hk_inp,k,H_real,H_compl)
!    !unfolds the real Hamiltonian into a 2 real Hamiltonians (one for the real and for the complex part)
!    type(t_h),intent(in)   :: Hk_inp
!    real(8),intent(in)          :: k(3)
!    class(t_h),allocatable,intent(out)     :: H_real,H_compl
!    class(t_h),allocatable                 :: Htmp_real,Htmp_compl
!
!    type(parameters_ham_init)   :: hinit
!    integer     :: i_ham,N_ham
!    complex(8),allocatable  ::  val(:)
!    integer,allocatable     ::  row(:),col(:)
!
!    real(8)  :: phase
!
!    Call set_H(H_real,h_io)
!    Call set_H(H_compl,h_io)
!    allocate(Htmp_real,mold=H)
!    allocate(Htmp_compl,mold=H)
!
!    N_ham=size(HK_inp%H)
!    do i_ham=1,N_ham
!        phase=dot_product(HK_inp%diffR(:,i_ham),k)
!        Call Hk_inp%H(i_ham)%get_hinit(hinit)
!        Call Hk_inp%H(i_ham)%get_par(val,row,col)
!
!        val=val*cmplx(cos(phase),sin(phase),8)
!
!        Call Htmp_real%init_coo(real(val,8),row,col,hinit)
!        Call Htmp_compl%init_coo(aimag(val,8),row,col,hinit)
!
!        Call H_real%add(Htmp_real)
!        Call H_compl%add(Htmp_compl)
!
!        Call Htmp_real%destroy()
!        Call Htmp_compl%destroy()
!        deallocate(val,row,col)
!    enddo
!end subroutine

end module m_construction_Hk
