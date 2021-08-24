module m_tb_k_zheevd
use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
#ifdef CPP_LAPACK

!use m_H_tb_public
use m_work_ham_single, only: work_ham, N_work
!use m_tb_k_base,  only: H_k_base
use m_tb_k_dense, only: H_k_dense
private
public :: H_k_zheevd

type,extends(H_k_dense) ::H_k_zheevd
    integer ::  size_work(3)
contains
    procedure   :: set_work     !sets the work arrays to the correct size
    procedure   :: calc_eval        !subroutine to calculate the eigenvalues using the internally set matrix
    procedure   :: calc_evec        !subroutine to calculate the eigenvectors and eigenvalues using the internally set matrix
end type

contains

subroutine set_work(this,work)
    class(H_k_zheevd),intent(inout)     :: this
    type(work_ham),intent(inout)        :: work
    integer                             :: sizes(N_work)

    complex(8)  :: cwork(1)
    real(8)     :: rwork(1)
    integer     :: iwork(1)
    integer     :: info

    sizes=0
    Call zheevd('V', 'U', this%dimH, this%H, this%dimH, this%eval,&
                & cwork, -1, rwork, -1, iwork, -1, info)
    sizes(1)=int(rwork(1))   !real
    sizes(2)=    iwork(1)    !int
    sizes(3)=int(cwork(1))   !complex
    this%size_work=sizes
    Call work%set(sizes)
end subroutine

subroutine calc_eval(this,Nin,eval,Nout,work)
    class(H_k_zheevd),intent(inout) :: this
    integer,intent(in)              :: Nin  !size of eigenvalue input array
    real(8),intent(inout)           :: eval(Nin)    !eigenvalue array
    integer,intent(out)             :: Nout !calculated number of eigenvalues
    type(work_ham),intent(inout)    :: work !work array that should be set to correct sizes
    !internal
    integer                 :: info
    external zheevd

    integer                 :: lcwork, lrwork, liwork

    lrwork=this%size_work(1)
    liwork=this%size_work(2)
    lcwork=this%size_work(3)

    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVALUE AS HAMILTONIAN IS NOT SET"
    if(Nin/=this%dimH) ERROR STOP "Size of eigenvalue array for tight-binding k-space Hamiltonian too small"
    Call zheevd('N', 'U', this%dimH, this%H, this%dimH, eval, &
                & work%cmplx_arr, lcwork, work%real_arr, lrwork, work%int_arr, liwork, info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    Nout=this%dimH
end subroutine

subroutine calc_evec(this,Nin,eval,evec,Nout,work)
    class(H_k_zheevd),intent(inout)      :: this
    integer,intent(in)                  :: Nin  !size of eigenvalue input array
    real(8),intent(inout)               :: eval(Nin)
    complex(8),intent(inout)            :: evec(this%dimH,Nin)
    integer,intent(out)                 :: Nout !calculated number of eigenvalues
    type(work_ham),intent(inout)        :: work !work array that should be set to correct sizes
    !internal
    integer                 :: info
    external zheevd

    integer                 :: lcwork, lrwork, liwork

    lrwork=this%size_work(1)
    liwork=this%size_work(2)
    lcwork=this%size_work(3)

    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVALUE AS HAMILTONIAN IS NOT SET"
    if(Nin/=this%dimH) ERROR STOP "Size of eigenvalue array for tight-binding k-space Hamiltonian too small"
    Call zheevd('V', 'U', this%dimH, this%H, this%dimH, eval, &
                & work%cmplx_arr, lcwork, work%real_arr, lrwork, work%int_arr, liwork, info)

    if(info/=0) ERROR STOP "LAPACK ERROR"
    Nout=this%dimH
    if(Nout>Nin)then
        write(error_unit,'(/A)') "Trying to diagonalize Hamiltonian with zheevd, but output Hamiltonian size has been set to low (TB_diag_estNe)"
        write(error_unit,'(AI8)') "Found eigenvalues in interval:", Nout
        write(error_unit,'(AI8)') "Assumend maximal array size  :", Nin
        ERROR STOP "INCREASE ARRAY SIZE (TB_diag_estNe) OR DECREASE ENERGY WINDOW TB_diag_Ebnd."
    endif
    evec=this%H
end subroutine

#endif
end module
