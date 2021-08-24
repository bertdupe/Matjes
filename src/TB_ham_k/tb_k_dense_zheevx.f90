module m_tb_k_zheevx
use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
#ifdef CPP_LAPACK

!use m_H_tb_public
use m_work_ham_single, only: work_ham, N_work
!use m_tb_k_base,  only: H_k_base
use m_tb_k_dense, only: H_k_dense
private
public :: H_k_zheevx

type,extends(H_k_dense) ::H_k_zheevx
    integer ::  size_work(3)
contains
    procedure   :: get_size_eval
    procedure   :: set_work    !sets the work arrays to the correct size
    procedure   :: calc_eval        !subroutine to calculate the eigenvalues using the internally set matrix
    procedure   :: calc_evec        !subroutine to calculate the eigenvectors and eigenvalues using the internally set matrix
end type

contains

integer function  get_size_eval(this)
    class(H_k_zheevx),intent(in)    :: this
    get_size_eval=this%Ne
end function

subroutine set_work(this,work)
    class(H_k_zheevx),intent(inout)     :: this
    type(work_ham),intent(inout)        :: work
    integer                             :: sizes(N_work)

    complex(8)  :: cwork(1)
    real(8)     :: rwork(1)
    integer     :: iwork(1)
    integer     :: Nout
    integer     :: info
    complex(8)  :: evec(1,1)
    integer     :: ifail(this%dimH)

    sizes=0
    Call zheevx('V', 'V', 'U', this%dimH, this%H, this%dimH, this%Ebnd(1), this%Ebnd(2),&
                & 0, 0, this%diag_acc, Nout, this%eval, evec, this%dimH,&
                & cwork, -1, rwork, iwork, ifail, info)
    sizes(1)=7*this%dimH    !real r-work
    sizes(2)=6*this%dimH    !integer, iwork and ifail
    sizes(3)=int(cwork(1))  !complex
    this%size_work=sizes
    Call work%set(sizes)
end subroutine

subroutine calc_eval(this,Nin,eval,Nout,work)
    class(H_k_zheevx),intent(inout) :: this
    integer,intent(in)              :: Nin  !size of eigenvalue input array
    real(8),intent(inout)           :: eval(Nin)    !eigenvalue array
    integer,intent(out)             :: Nout !calculated number of eigenvalues
    type(work_ham),intent(inout)    :: work !work array that should be set to correct sizes
    !internal
    integer                 :: info
    complex(8)              :: evec(1,1)
    external zheevx

    integer                 :: lcwork

    lcwork=this%size_work(3)

    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVALUE AS HAMILTONIAN IS NOT SET"
    Call zheevx('N', 'V', 'U', this%dimH, this%H, this%dimH, this%Ebnd(1), this%Ebnd(2),&
                & 0, 0, this%diag_acc, Nout, eval, evec, this%dimH,&
                & work%cmplx_arr, lcwork, work%real_arr, work%int_arr(1:this%dimH*5), work%int_arr(this%dimH*5+1:this%dimH*6), info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    if(Nout>Nin)then
        write(error_unit,'(/A)') "Trying to diagonalize Hamiltonian with zheevx, but output Hamiltonian size has been set to low (TB_diag_estNe)"
        write(error_unit,'(AI8)') "Found eigenvalues in interval:", Nout
        write(error_unit,'(AI8)') "Assumend maximal array size  :", Nin
        ERROR STOP "INCREASE ARRAY SIZE (TB_diag_estNe) OR DECREASE ENERGY WINDOW TB_diag_Ebnd."
    endif
end subroutine

subroutine calc_evec(this,Nin,eval,evec,Nout,work)
    class(H_k_zheevx),intent(inout)      :: this
    integer,intent(in)                  :: Nin  !size of eigenvalue input array
    real(8),intent(inout)               :: eval(Nin)
    complex(8),intent(inout)            :: evec(this%dimH,Nin)
    integer,intent(out)                 :: Nout !calculated number of eigenvalues
    type(work_ham),intent(inout)        :: work !work array that should be set to correct sizes
    !internal
    integer                 :: info
    external zheevx

    integer                 :: lcwork

    lcwork=this%size_work(3)

    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVALUE AS HAMILTONIAN IS NOT SET"
    Call zheevx('V', 'V', 'U', this%dimH, this%H, this%dimH, this%Ebnd(1), this%Ebnd(2),&
                & 0, 0, this%diag_acc, Nout, eval, evec, this%dimH,&
                & work%cmplx_arr, lcwork, work%real_arr, work%int_arr(1:this%dimH*5), work%int_arr(this%dimH*5+1:this%dimH*6), info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    if(Nout>Nin)then
        write(error_unit,'(/A)') "Trying to diagonalize Hamiltonian with zheevx, but output Hamiltonian size has been set to low (TB_diag_estNe)"
        write(error_unit,'(AI8)') "Found eigenvalues in interval:", Nout
        write(error_unit,'(AI8)') "Assumend maximal array size  :", Nin
        ERROR STOP "INCREASE ARRAY SIZE (TB_diag_estNe) OR DECREASE ENERGY WINDOW TB_diag_Ebnd."
    endif
end subroutine

#endif
end module
