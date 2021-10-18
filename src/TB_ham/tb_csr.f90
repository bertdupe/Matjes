module m_H_tb_csr
#ifdef CPP_MKL
use m_H_tb_base
use m_H_tb_coo
private
public H_feast_csr

type,extends(H_TB_coo_based),abstract :: H_tb_csr
    private
    integer, allocatable    :: i(:),j(:) 
    complex(8), allocatable :: csr(:)
    contains
    procedure   :: set_from_Hcoo
    procedure   :: add_child
    procedure   :: copy_child
    procedure   :: destroy_child
    procedure   :: mv
    procedure   :: mult_r
end type

type,extends(H_TB_csr)  :: H_feast_csr
    contains
    procedure   :: get_evec => evec_feast
end type


contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  GENERAL ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_from_Hcoo(this,H_coo)
    class(H_tb_csr),intent(inout)   :: this
    type(H_tb_coo),intent(inout)    :: H_coo
    integer             :: nnz
    complex(8),allocatable :: val_coo(:)
    integer,allocatable :: rowind(:),colind(:)
    integer,parameter   ::  job(6)=[2,1,1,0,0,0]
    integer             ::  info

    Call this%init_otherH(H_coo)
    Call H_coo%pop_par(nnz,val_coo,rowind,colind)
    allocate(this%csr(nnz),this%j(nnz),this%i(this%dimH+1))
    call mkl_zcsrcoo(job, this%dimH, this%csr, this%j, this%i, nnz, val_coo, rowind, colind, info)
    if(info/=0) STOP "mkl_zcsrcoo error"
end subroutine 

subroutine mv(this,Hout)
    class(H_tb_csr),intent(inout) :: this
    class(H_TB),intent(inout)       :: Hout
    
    select type(Hout)
    class is(H_TB_csr)
        Call Hout%init_otherH(this)
        Call move_alloc(this%i,Hout%i)
        Call move_alloc(this%j,Hout%j)
        Call move_alloc(this%csr,Hout%csr)
    class default
        STOP "Cannot move H_TB_csr type to Hamiltonian that is not a class of H_TB_csr"
    end select
    Call this%destroy()
end subroutine

subroutine copy_child(this,Hout)
    class(H_tb_csr),intent(in)  :: this
    class(H_TB),intent(inout)   :: Hout
    
    select type(Hout)
    class is(H_tb_csr)
        allocate(Hout%csr,source=this%csr)
        allocate(Hout%i,source=this%i)
        allocate(Hout%j,source=this%j)
    class default
        STOP "Cannot copy H_tb_csr type with Hamiltonian that is not a class of H_tb_csr"
    end select
end subroutine

recursive subroutine add_child(this,H_in)
    class(H_tb_csr),intent(inout)   :: this
    class(H_tb),intent(in)          :: H_in

    class(H_tb_csr),allocatable     :: Htmp
    type(H_tb_coo)                  :: Hcoo_tmp

    integer, allocatable    :: ia(:),ja(:) 
    complex(8), allocatable :: acsr(:)
    complex(8),parameter    :: beta=(1.0d0,0.0d0)
    integer                 :: info
    complex(8)              :: cdummy(1)
    integer                 :: idummy(1)

    select type(H_in)
    class is(H_tb_csr)
        Call move_alloc(this%i,ia)
        Call move_alloc(this%j,ja)
        Call move_alloc(this%csr,acsr)
        allocate(this%i(this%dimH+1))
        call mkl_zcsradd('n', 1, 0,this%dimH ,this%dimH, acsr ,ja, ia, beta, H_in%csr, H_in%j, H_in%i, cdummy  , idummy, this%i, 0, info)
        if(info/=0) STOP "ERROR in request 1 of mkl_zcsradd"
        allocate(this%j  (this%i(this%dimH+1)-1))
        allocate(this%csr(this%i(this%dimH+1)-1))
        call mkl_zcsradd('n', 2, 0,this%dimH ,this%dimH, acsr ,ja, ia, beta, H_in%csr, H_in%j, H_in%i, this%csr, this%j, this%i, 0, info)
        if(info/=0) STOP "ERROR in request 2 of mkl_zcsradd"
    class is(H_TB_coo)
        Call H_in%copy(Hcoo_tmp)    !make local copy since add should not destroy the added array
        allocate(Htmp,mold=this)    
        Call Htmp%set_from_Hcoo(Hcoo_tmp)
        Call this%add(Htmp)
        Call Htmp%destroy()
        Call Hcoo_tmp%destroy()
    class default
        ERROR STOP "Cannot add H_tb_csr type with Hamiltonian that is not a class of H_tb_csr"
    end select
end subroutine

subroutine destroy_child(this)
    class(H_tb_csr),intent(inout)    :: this

    if(allocated(this%csr)) deallocate(this%csr)
    if(allocated(this%i))   deallocate(this%i  )
    if(allocated(this%j))   deallocate(this%j  )
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  multiplication routine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mult_r(this,vec,res,alpha,beta)
    class(H_tb_csr),intent(in)   :: this
    complex(8),intent(in   )        :: vec(this%dimH)
    complex(8),intent(inout)        :: res(this%dimH)
    complex(8),intent(in),optional  :: alpha
    complex(8),intent(in),optional  :: beta
    ERROR STOP "IMPLEMENT MULT_R TB_CSR"
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  FEAST ROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine evec_feast(this,eval,evec)
    class(H_feast_csr),intent(in)       ::  this
    real(8),intent(out),allocatable     ::  eval(:)
    complex(8),intent(out),allocatable  ::  evec(:,:)

    integer                     :: fpm(128)
    real(8),allocatable         :: e(:)
    complex(8),allocatable      :: x(:,:)
    real(8)                     :: emin,emax
    real(8)                     :: epsout
    integer                     :: loop
    integer                     :: m0,m
    real(8),allocatable         :: res(:)
    integer                     :: info

    Call feastinit(fpm) 
    fpm(1)=1
    fpm(2)=8
    fpm(3)=-nint(log10(this%diag_acc))
    emin=this%Ebnd(1)
    emax=this%Ebnd(2)
    m0=this%estNe
    if(m0==0.or.m0>this%dimH) m0=this%dimH
    allocate(e(m0),source=0.0d0)
    allocate(x(this%dimh,m0),source=(0.0d0,0.0d0))
    allocate(res(this%dimH),source=0.0d0)
    call zfeast_hcsrev('F', this%dimH, this%csr, this%i, this%j, fpm, epsout, loop, emin, emax, m0, e, x, m, res, info)
    if(info/=0) STOP "sparse diagonalization failed"
    allocate(eval,source=e(1:m))
    deallocate(e)
    allocate(evec,source=x(1:this%dimH,1:m))
    deallocate(x,res)
end subroutine
#endif
end module
