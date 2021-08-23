module m_tb_k_dense
use m_TB_types, only: parameters_TB_IO_H
use m_Hk, only: Hk_inp_t
use m_H_tb_public
use m_work_ham_single, only: work_ham, N_work
use m_tb_k_base, only: H_k_base

private
public H_k_dense
type,extends(H_k_base),abstract  ::H_k_dense 
    complex(8),allocatable          :: H(:,:)       !temporary Hamiltonian array which is set and evaluated [dimH,dimH]         !put into work array?
    real(8),allocatable             :: eval(:)      !temporary array for eigenvalues                                            !put into work array?
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
    type(H_TB_coo)     :: H_tmp
#ifdef CPP_LAPACK
    type(H_zheev)      :: H_loc    !real space Hamiltonian as the given distance by diffR [N_H]

    integer     ::  iH

    Call this%init_base(Hk_inp,io)

    !get real-space Hamiltonians (H_R) from Hcoo input
    allocate(this%H_R(this%dimH,this%dimH,this%N_H))
    do iH=1,size(Hk_inp%H)
        Call Hk_inp%H(iH)%copy(H_tmp)
        Call H_loc%set_from_Hcoo(H_tmp)
        Call H_loc%get_H(this%H_R(:,:,iH))
        Call H_tmp%destroy()
        Call H_loc%destroy()
    enddo
#else
    ERROR STOP "CANNOT INITIALIZE H_K_dense without CPP_LAPACK"
#endif

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
    integer ::  dimH2
    dimH2=this%dimH*this%dimH
    Call get_Hk(dimH2,this%N_H,k,this%H,this%H_r,this%diffR)
#endif
end subroutine

subroutine get_Hk(dimH2,Nr,k,H,H_r,diffR)
    integer,intent(in)              :: dimH2, Nr
    real(8),intent(in)              :: k(3)
    real(8),intent(in)              :: diffR(3,Nr)
    complex(8),intent(out)          :: H(dimH2)
    complex(8),intent(in)           :: H_r(dimH2,Nr)
    real(8)     :: phase_r(Nr)
    complex(8)  :: phase_c(Nr)
    integer     :: iH
#ifdef CPP_BLAS
    external zgemv
#endif

    phase_r=matmul(k,diffR)
    phase_c=cmplx(cos(phase_r),sin(phase_r),8)
#ifdef CPP_BLAS
    H=(0d0,0d0)
    Call zgemv('n', dimH2, Nr, (1.0d0,0.0d0), H_r,dimH2, phase_c ,1 ,(0.0d0,0d0),H,1)
#else
    H=matmul(H_r,phase_c)
#endif
end subroutine
end module
