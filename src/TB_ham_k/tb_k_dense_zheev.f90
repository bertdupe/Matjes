module m_tb_k_zheev
#ifdef CPP_LAPACK

!use m_H_tb_public
use m_work_ham_single, only: work_ham, N_work
!use m_tb_k_base,  only: H_k_base
use m_tb_k_dense, only: H_k_dense
private
public :: H_k_zheev

type,extends(H_k_dense) ::H_k_zheev
contains
    procedure   ::  set_work    !sets the work arrays to the correct size
    procedure   ::  calc_eval        !subroutine to calculate the eigenvalues using the internally set matrix
    procedure   ::  calc_evec        !subroutine to calculate the eigenvectors and eigenvalues using the internally set matrix
end type

contains

subroutine set_work(this,work)
    class(H_k_zheev),intent(inout)      :: this
    type(work_ham),intent(inout)        :: work
    integer                             :: sizes(N_work)

    sizes=0
    sizes(1)=max(1,3*this%dimH-2)   !real(RWORK)
    sizes(3)=max(1,2*this%dimH-1)   !complex(WORK)
    Call work%set(sizes)
end subroutine

subroutine calc_eval(this,Nin,eval,Nout,work)
    class(H_k_zheev),intent(inout)  :: this
    integer,intent(in)              :: Nin  !size of eigenvalue input array
    real(8),intent(inout)           :: eval(Nin)    !eigenvalue array
    integer,intent(out)             :: Nout !calculated number of eigenvalues
    type(work_ham)                  :: work !work array that should be set to correct sizes
    !internal
    integer                 :: info,lcwork,lrwork
    external zheev

    lcwork=max(1,2*this%dimH-1)
    lrwork=max(1,3*this%dimH-2)
#ifdef CPP_DEBUG
    if(size(work%real_arr )/=lrwork) ERROR STOP "REAL WORK SIZE SET INCORRECTLY"
    if(size(work%cmplx_arr)/=lcwork) ERROR STOP "COMPLEX WORK SIZE SET INCORRECTLY"
#endif
    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVALUE AS HAMILTONIAN IS NOT SET"
    Call zheev('N', 'U', this%dimH, this%H, this%dimH, this%eval, work%cmplx_arr(1:lcwork), lcwork, work%real_arr(1:lrwork), info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    if(Nin<this%dimH)then
        eval=this%eval(1:Nin)
        Nout=Nin
    else
        Nout=this%dimH
        eval(1:Nout)=this%eval
    endif
end subroutine

subroutine calc_evec(this,Nin,eval,evec,Nout,work)
    class(H_k_zheev),intent(inout)      :: this
    integer,intent(in)                  :: Nin  !size of eigenvalue input array
    real(8),intent(inout)               :: eval(Nin)
    complex(8),intent(inout)            :: evec(this%dimH,Nin)
    integer,intent(out)                 :: Nout !calculated number of eigenvalues
    type(work_ham)                      :: work !work array that should be set to correct sizes
    !internal
    integer                 :: info,lcwork,lrwork
    external zheev
    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVALUE AS HAMILTONIAN IS NOT SET"
    lcwork=max(1,2*this%dimH-1)
    lrwork=max(1,3*this%dimH-2)
#ifdef CPP_DEBUG
    if(size(work%real_arr )/=lrwork) ERROR STOP "REAL WORK SIZE SET INCORRECTLY"
    if(size(work%cmplx_arr)/=lcwork) ERROR STOP "COMPLEX WORK SIZE SET INCORRECTLY"
#endif
    Call zheev('V', 'U', this%dimH, this%H, this%dimH, this%eval, work%cmplx_arr(1:lcwork), lcwork, work%real_arr(1:lrwork), info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    if(Nin<this%dimH)then
        Nout=Nin
        eval=this%eval(1:Nout)
        evec=this%H(1:this%dimH,1:Nout)
    else
        Nout=this%dimH
        eval(1:Nout)=this%eval
        evec(1:this%dimH,1:Nout)=this%H
    endif
end subroutine

#endif
end module
