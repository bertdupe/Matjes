module m_tb_k_feast
use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
#ifdef CPP_MKL

!use m_H_tb_public
use m_work_ham_single, only: work_ham, N_work
!use m_tb_k_base,  only: H_k_base
use m_tb_k_dense, only: H_k_dense
private
public :: H_k_feast

type,extends(H_k_dense) ::H_k_feast
    integer :: size_work(3)
    integer :: fpm(128)
contains
    procedure   :: get_size_eval
    procedure   :: set_work    !sets the work arrays to the correct size
    procedure   :: calc_eval        !subroutine to calculate the eigenvalues using the internally set matrix
    procedure   :: calc_evec        !subroutine to calculate the eigenvectors and eigenvalues using the internally set matrix
end type

contains

integer function  get_size_eval(this)
    class(H_k_feast),intent(in)    :: this
    get_size_eval=this%Ne
end function

subroutine set_work(this,work)
    class(H_k_feast),intent(inout)     :: this
    type(work_ham),intent(inout)        :: work
    integer                             :: sizes(N_work)

    Call feastinit(this%fpm) 
    this%fpm(1)=0
    this%fpm(2)=8
    this%fpm(3)=-nint(log10(this%diag_acc))

    sizes=0
    sizes(1)=this%Ne
    this%size_work=sizes
    Call work%set(sizes)
end subroutine

subroutine calc_eval(this,Nin,eval,Nout,work)
    class(H_k_feast),intent(inout) :: this
    integer,intent(in)              :: Nin  !size of eigenvalue input array
    real(8),intent(inout)           :: eval(Nin)    !eigenvalue array
    integer,intent(out)             :: Nout !calculated number of eigenvalues
    type(work_ham),intent(inout)    :: work !work array that should be set to correct sizes
    !internal

    real(8)                     :: epsout
    integer                     :: loop
    integer                     :: info
    complex(8),allocatable      :: x(:,:)

    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVALUE AS HAMILTONIAN IS NOT SET"
    allocate(x(this%dimH,Nin))
    call zfeast_heev ( 'F', this%dimH, this%H, this%dimH, this%fpm, epsout, loop, this%Ebnd(1), this%Ebnd(2), Nin, eval, x, Nout, work%real_arr, info )
    if(info/=0) ERROR STOP "Feast ERROR"
    if(Nout>Nin)then
        write(error_unit,'(/A)') "Trying to diagonalize Hamiltonian with feast, but output Hamiltonian size has been set to low (TB_diag_estNe)"
        write(error_unit,'(AI8)') "Found eigenvalues in interval:", Nout
        write(error_unit,'(AI8)') "Assumend maximal array size  :", Nin
        ERROR STOP "INCREASE ARRAY SIZE (TB_diag_estNe) OR DECREASE ENERGY WINDOW TB_diag_Ebnd."
    endif
end subroutine

subroutine calc_evec(this,Nin,eval,evec,Nout,work)
    class(H_k_feast),intent(inout)      :: this
    integer,intent(in)                  :: Nin  !size of eigenvalue input array
    real(8),intent(inout)               :: eval(Nin)
    complex(8),intent(inout)            :: evec(this%dimH,Nin)
    integer,intent(out)                 :: Nout !calculated number of eigenvalues
    type(work_ham),intent(inout)        :: work !work array that should be set to correct sizes
    !internal
    real(8)                     :: epsout
    integer                     :: loop
    integer                     :: info

    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVALUE AS HAMILTONIAN IS NOT SET"
    call zfeast_heev ( 'F', this%dimH, this%H, this%dimH, this%fpm, epsout, loop, this%Ebnd(1), this%Ebnd(2), Nin, eval, evec, Nout, work%real_arr, info )
    if(info/=0) ERROR STOP "Feast ERROR"
    if(Nout>Nin)then
        write(error_unit,'(/A)') "Trying to diagonalize Hamiltonian with feast, but output Hamiltonian size has been set to low (TB_diag_estNe)"
        write(error_unit,'(AI8)') "Found eigenvalues in interval:", Nout
        write(error_unit,'(AI8)') "Assumend maximal array size  :", Nin
        ERROR STOP "INCREASE ARRAY SIZE (TB_diag_estNe) OR DECREASE ENERGY WINDOW TB_diag_Ebnd."
    endif
end subroutine

#endif
end module
