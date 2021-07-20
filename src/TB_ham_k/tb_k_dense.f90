module m_tb_k_dense
use m_TB_types, only: parameters_TB_IO_H
use m_Hk, only: Hk_inp_t
use m_H_tb_public
use m_work_ham_single, only: work_ham, N_work
use m_tb_k_base, only: H_k_base

private
public H_k_dense
type,extends(H_k_base),abstract  ::H_k_dense 
    complex(8),allocatable          :: H(:,:)       !temporary Hamiltonian array which is set and evaluated [dimH,dimH]
    real(8),allocatable             :: eval(:)      !temporary array for eigenvalues
    complex(8),allocatable          :: H_R(:,:,:)   !real-space Hamiltonians as the given distance by diffR [dimH,dimH,N_H] 
contains
    procedure           :: init
    procedure           :: set_k
end type
contains
subroutine init(this,Hk_inp,io)
    class(H_k_dense),intent(inout)      :: this
    type(Hk_inp_t),intent(inout)        :: Hk_inp
    type(parameters_TB_IO_H),intent(in) :: io
    type(H_TB_coo)                      :: H_tmp
    type(H_zheev)                       :: H_loc    !real space Hamiltonian as the given distance by diffR [N_H]

    integer     ::  iH

    allocate(this%diffR,source=Hk_inp%diffR)
    this%dimH=Hk_inp%H(1)%dimH
    this%N_H=size(this%diffR,2)

    !get real-space Hamiltonians (H_R) from Hcoo input
    allocate(this%H_R(this%dimH,this%dimH,this%N_H))
    do iH=1,size(Hk_inp%H)
        Call Hk_inp%H(iH)%copy(H_tmp)
        Call H_loc%set_from_Hcoo(H_tmp)
        Call H_loc%get_H(this%H_R(:,:,iH))
        Call H_tmp%destroy()
        Call H_loc%destroy()
    enddo

    !prepare temp H
    allocate(this%H(this%dimH,this%dimH))
    allocate(this%eval(this%dimH))

    Call this%do_set(.true.)
end subroutine

subroutine set_k(this,k) 
    class(H_k_dense),intent(inout)      :: this
    real(8),intent(in)                  :: k(3)

#if 0
    real(8)     :: phase_r
    complex(8)  :: phase_c
    integer  :: iH
    this%H=(0.0d0,0.0d0)
    do iH=1,this%N_H
        phase_r=dot_product(this%diffR(:,iH),k)
        phase_c=cmplx(cos(phase_r),sin(phase_r),8)
        this%H=this%H+phase_c*this%H_R(:,:,iH)
    enddo
#else
    Call get_Hk(this%dimH,this%N_H,k,this%H,this%H_r,this%diffR)
#endif
end subroutine

subroutine get_Hk(dimH,Nr,k,H,H_r,diffR)
    integer,intent(in)              :: dimH, Nr
    real(8),intent(in)              :: k(3)
    real(8),intent(in)              :: diffR(3,Nr)
    complex(8),intent(out)          :: H(dimH*dimH)
    complex(8),intent(in)           :: H_r(dimH*dimH,Nr)
    real(8)     :: phase_r(Nr)
    complex(8)  :: phase_c(Nr)
    integer     :: iH

    phase_r=matmul(k,diffR)
    phase_c=cmplx(cos(phase_r),sin(phase_r),8)
    H=matmul(H_r,phase_c)
end subroutine

!#ifdef CPP_LAPACK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!  ZHEEV ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!subroutine set_work_zheev(this,work)
!    class(H_k_dense),intent(inout)      :: this
!    type(work_ham),intent(inout)        :: work
!    integer                             :: sizes(N_work)
!
!    sizes=0
!    sizes(1)=max(1,3*this%dimH-2)   !real(RWORK)
!    sizes(3)=max(1,2*this%dimH-1)   !complex(WORK)
!    Call work%set(sizes)
!end subroutine
!
!subroutine eval_zheev(this,Nin,eval,Nout,work)
!    class(H_zheev),intent(inout)    :: this
!    integer,intent(in)              :: Nin  !size of eigenvalue input array
!    real(8),intent(inout)           :: eval(Nin)    !eigenvalue array
!    integer,intent(out)             :: Nout !calculated number of eigenvalues
!    type(work_ham)                  :: work !work array that should be set to correct sizes
!    !internal
!    integer                 :: info,lcwork,lrwork
!    external zheev
!
!    lcwork=max(1,2*this%dimH-1)
!    rwork=max(1,3*this%dimH-2)
!#ifdef CPP_DEBUG
!    if(size(work%real_arr )/=lrwork) ERROR STOP "REAL WORK SIZE SET INCORRECTLY"
!    if(size(work%cmplx_arr)/=lwork) ERROR STOP "COMPLEX WORK SIZE SET INCORRECTLY"
!#endif
!    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVALUE AS HAMILTONIAN IS NOT SET"
!    Call zheev('N', 'U', this%dimH, this%H, this%dimH, this%eval, work%cmplx_arr(1:lcwork), lcwork, work%real_arr(1:lrwork), info)
!    if(info/=0) ERROR STOP "LAPACK ERROR"
!    if(Nin<this%dimH)then
!        eval=this%eval(1:Nin)
!        Nout=Nin
!    else
!        Nout=this%dimH
!        eval(1:Nout)=this%eval
!    endif
!end subroutine
!
!subroutine evec_zheev(this,Nin,eval,evec,Nout,work)
!    class(H_zheev),intent(inout)        :: this
!    integer,intent(in)                  :: Nin  !size of eigenvalue input array
!    real(8),intent(inout)               :: eval(Nin)
!    complex(8),intent(inout)            :: evec(this%dimH,Nin)
!    integer,intent(out)                 :: Nout !calculated number of eigenvalues
!    type(work_ham)                      :: work !work array that should be set to correct sizes
!    !internal
!    integer                 :: info,lcwork,lrwork
!    external zheev
!    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVALUE AS HAMILTONIAN IS NOT SET"
!    lcwork=max(1,2*this%dimH-1)
!    rwork=max(1,3*this%dimH-2)
!#ifdef CPP_DEBUG
!    if(size(work%real_arr )/=lrwork) ERROR STOP "REAL WORK SIZE SET INCORRECTLY"
!    if(size(work%cmplx_arr)/=lwork) ERROR STOP "COMPLEX WORK SIZE SET INCORRECTLY"
!#endif
!    Call zheev('V', 'U', this%dimH, this%H, this%dimH, this%eval, work%cmplx_arr(1:lcwork), lcwork, work%real_arr(1:lrwork), info)
!    if(info/=0) ERROR STOP "LAPACK ERROR"
!    if(Nin<this%dimH)then
!        Nout=Nin
!        eval=this%eval(1:Nout)
!        evec=this%H(1:this%dimH,1:Nout)
!    else
!        Nout=this%dimH
!        eval(1:Nout)=this%eval
!        evec(1:this%dimH,1:Nout)=this%H
!    endif
!    evec=this%H
!end subroutine
!#endif

end module
