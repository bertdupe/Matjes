module m_FT_Ham_dense_base
use m_H_type
use m_parameters_FT_Ham
use m_io_files_utils
use m_FT_Ham_base
use m_FT_Ham_coo
!use m_FT_Ham_zheev
use m_FT_Ham_dense
implicit none

private
public FT_H_dense_base

type,extends(t_H_base),abstract     :: FT_H_dense_base
    complex(8),allocatable          :: eval(:)      !temporary array for eigenvalues
contains
    procedure           :: init
    procedure           :: set_k
end type

contains

subroutine init(this,Hk_inp)
    class(FT_H_dense_base),intent(inout)      :: this
    type(H_inp_real_to_k),intent(inout)       :: Hk_inp
    type(H_inp_k_coo)                         :: H_tmp
!    type(FT_Ham_zheev)                        :: H_loc    !real space Hamiltonian as the given distance by diffR [N_H]

    integer     ::  iH,io_input

    allocate(this%diffR,source=Hk_inp%diffR)
    this%dimH=Hk_inp%H(1)%dimH
    this%N_H=size(this%diffR,2)

    !get real-space Hamiltonians (H_R) from Hcoo input
    allocate(this%H_R(this%dimH,this%dimH,this%N_H))
    do iH=1,size(Hk_inp%H)
        Call Hk_inp%H(iH)%copy(H_tmp)
!        Call H_loc%set_from_Hcoo(H_tmp)
!        Call H_loc%get_H(this%H_R(:,:,iH))
        Call H_tmp%destroy()
!        Call H_loc%destroy()
    enddo

    ! initialize diagonalization parameters
    if (this%dimH(1).eq.this%dimH(2)) this%io_H%dimH=this%dimH(1)
    io_input=open_file_read('input')
    call this%io_H%read_file(this,io_input,'input')
    call close_file('input',io_input)

    !prepare temp H
    allocate(this%H(this%dimH,this%dimH))
    allocate(this%eval(this%dimH))

    Call this%do_set(.true.)
end subroutine

subroutine set_k(this,k)
    class(FT_H_dense_base),intent(inout)  :: this
    real(8),intent(in)                    :: k(3)

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
end module
