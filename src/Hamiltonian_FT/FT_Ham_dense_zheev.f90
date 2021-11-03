module m_FT_Ham_zheev
#ifdef CPP_LAPACK
use m_work_ham_single, only: work_ham, N_work
use m_FT_Ham_dense
use m_parameters_FT_Ham
implicit none

private
public :: FT_Ham_zheev

type,extends(FT_H_dense) :: FT_Ham_zheev
    complex(8), allocatable    :: cmplx_arr(:)
    real(8), allocatable       :: real_arr(:)
contains
    procedure   ::  set_work         !sets the work arrays to the correct size
    procedure   ::  calc_eval        !subroutine to calculate the eigenvalues using the internally set matrix
    procedure   ::  calc_evec        !subroutine to calculate the eigenvectors and eigenvalues using the internally set matrix
end type

contains

subroutine set_work(this,eval,evec)
    class(FT_Ham_zheev),intent(inout)     :: this
    real(8),allocatable,intent(out)       :: eval(:)      ! array containing the eigenvalues
    complex(8),allocatable,intent(out)    :: evec(:,:)   ! array containing the eigenvectors

    integer                             :: n_real,n_cmplx

    n_real=max(1,3*this%io_H%dimH-2)   !real(RWORK)
    n_cmplx=max(1,2*this%io_H%dimH-1)   !complex(WORK)

    allocate(eval(this%io_H%dimH),source=0.0d0)
    allocate(evec(this%io_H%dimH,this%io_H%dimH),source=(0.0d0,0.0d0))

    if(n_real>0) allocate(this%real_arr(n_real),source=0.0d0)
    if(n_cmplx>0) allocate(this%cmplx_arr(n_cmplx),source=(0.0d0,0.0d0))

end subroutine

subroutine calc_eval(this,Nin,eval,Nout)
    class(FT_Ham_zheev),intent(inout)  :: this
    integer,intent(in)                 :: Nin  !size of eigenvalue input array
    real(8),intent(inout)              :: eval(Nin)    !eigenvalue array
    integer,intent(out)                :: Nout !calculated number of eigenvalues
    !internal
    integer                 :: info,lcwork,lrwork
    complex(8), allocatable :: eval_tmp(:)
    external zheev

    allocate(eval_tmp(this%io_H%dimH),source=(0.0d0,0.0d0))
    lcwork=max(1,2*this%io_H%dimH-1)
    lrwork=max(1,3*this%io_H%dimH-2)
#ifdef CPP_DEBUG
    if(size(work%real_arr )/=lrwork) ERROR STOP "REAL WORK SIZE SET INCORRECTLY"
    if(size(work%cmplx_arr)/=lcwork) ERROR STOP "COMPLEX WORK SIZE SET INCORRECTLY"
#endif
    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVALUE AS HAMILTONIAN IS NOT SET"
    Call zheev('N', 'U', this%io_H%dimH, this%Hk, this%io_H%dimH, eval_tmp, this%cmplx_arr(1:lcwork), lcwork, this%real_arr(1:lrwork), info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    if(Nin<this%io_H%dimH)then
        eval=real(eval_tmp(1:Nin),8)
        Nout=Nin
    else
        Nout=this%io_H%dimH
        eval(1:Nout)=real(eval_tmp,8)
    endif
end subroutine

subroutine calc_evec(this,Nin,eval,evec,Nout)
    class(FT_Ham_zheev),intent(inout)   :: this
    integer,intent(in)                  :: Nin  !size of eigenvalue input array
    complex(8),intent(inout)            :: eval(Nin)
    complex(8),intent(inout)            :: evec(this%io_H%dimH,Nin)
    integer,intent(out)                 :: Nout !calculated number of eigenvalues
    !internal
    integer                 :: info,lcwork,lrwork
    external zheev

    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVALUE AS HAMILTONIAN IS NOT SET"
    lcwork=max(1,2*this%io_H%dimH-1)
    lrwork=max(1,3*this%io_H%dimH-2)
#ifdef CPP_DEBUG
    if(size(this%real_arr )/=lrwork) ERROR STOP "REAL WORK SIZE SET INCORRECTLY"
    if(size(this%cmplx_arr)/=lcwork) ERROR STOP "COMPLEX WORK SIZE SET INCORRECTLY"
#endif
    Call zheev('V', 'U', this%io_H%dimH, this%Hk, this%io_H%dimH, eval, this%cmplx_arr(1:lcwork), lcwork, this%real_arr(1:lrwork), info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    if(Nin<this%io_H%dimH)then
        Nout=Nin
        eval=eval(1:Nout)
        evec=this%Hk(1:this%io_H%dimH,1:Nout)
    else
        Nout=this%io_H%dimH
        eval(1:Nout)=eval
        evec(1:this%io_H%dimH,1:Nout)=this%Hk
    endif
end subroutine

#endif
end module
