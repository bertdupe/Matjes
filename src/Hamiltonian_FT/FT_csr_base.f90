module m_FT_csr
#ifdef CPP_MKL
use m_FT_Ham_coo
use m_FT_Ham_base
use m_H_combined
use m_parameters_FT_Ham
use m_H_type
private
public FT_H_feast_csr

type,extends(H_inp_k_coo) :: FT_H_csr
    private
    type(parameters_FT_HAM_IO)  :: ham_io
    integer, allocatable    :: i(:),j(:)
    complex(8), allocatable :: csr(:)
    contains
    procedure   :: set_from_Hcoo
    procedure   :: mv
    procedure   :: add_coo
end type

type,extends(FT_H_csr) :: FT_H_feast_csr
    contains
    procedure   :: get_evec => evec_feast
end type


contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  GENERAL ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_from_Hcoo(this,H_coo)
    class(FT_H_csr),intent(inout)       :: this
    class(H_inp_k_coo),intent(inout)    :: H_coo
    integer             :: nnz
    complex(8),allocatable :: val_coo(:)
    integer,allocatable :: rowind(:),colind(:)
    integer,parameter   ::  job(6)=[2,1,1,0,0,0]
    integer             ::  info

    Call this%init_otherH(H_coo)
    Call H_coo%pop_par(nnz,val_coo,rowind,colind)
    allocate(this%csr(nnz),this%j(nnz),this%i(this%ham_io%dimH+1))
    call mkl_zcsrcoo(job, this%ham_io%dimH, this%csr, this%j, this%i, nnz, val_coo, rowind, colind, info)
    if(info/=0) STOP "mkl_zcsrcoo error"
end subroutine

subroutine mv(this,Hout)
    class(FT_H_csr),intent(inout) :: this
    class(t_H_base),intent(inout)       :: Hout

    select type(Hout)
    class is(FT_H_csr)
        Call Hout%init_otherH(this)
        Call move_alloc(this%i,Hout%i)
        Call move_alloc(this%j,Hout%j)
        Call move_alloc(this%csr,Hout%csr)
    class default
        STOP "Cannot move FT_H_csr type to Hamiltonian that is not a class of H_TB_csr"
    end select
    Call this%destroy()
end subroutine

recursive subroutine add_coo(this,H_in)
    class(FT_H_csr),intent(inout)   :: this
    class(t_H_base),intent(in)      :: H_in

    class(FT_H_csr),allocatable     :: Htmp
    class(H_inp_k_coo),allocatable  :: Hcoo_tmp

    integer, allocatable    :: ia(:),ja(:)
    complex(8), allocatable :: acsr(:)
    complex(8),parameter    :: beta=(1.0d0,0.0d0)
    integer                 :: info
    complex(8)              :: cdummy(1)
    integer                 :: idummy(1)

    select type(H_in)
    class is(FT_H_csr)
        Call move_alloc(this%i,ia)
        Call move_alloc(this%j,ja)
        Call move_alloc(this%csr,acsr)
        allocate(this%i(this%ham_io%dimH+1))
        call mkl_zcsradd('n', 1, 0,this%ham_io%dimH ,this%ham_io%dimH, acsr ,ja, ia, beta, H_in%csr, H_in%j, H_in%i, cdummy  , idummy, this%i, 0, info)
        if(info/=0) STOP "ERROR in request 1 of mkl_zcsradd"
        allocate(this%j  (this%i(this%ham_io%dimH+1)-1))
        allocate(this%csr(this%i(this%ham_io%dimH+1)-1))
        call mkl_zcsradd('n', 2, 0,this%ham_io%dimH ,this%ham_io%dimH, acsr ,ja, ia, beta, H_in%csr, H_in%j, H_in%i, this%csr, this%j, this%i, 0, info)
        if(info/=0) STOP "ERROR in request 2 of mkl_zcsradd"
    class is(H_inp_k_coo)
        Call H_in%copy(Hcoo_tmp)    !make local copy since add should not destroy the added array
        allocate(Htmp,mold=this)
        Call Htmp%set_from_Hcoo(Hcoo_tmp)
        Call this%add_coo(Htmp)
        Call Htmp%destroy()
        Call Hcoo_tmp%destroy()
    class default
        ERROR STOP "Cannot add FT_H_csr type with Hamiltonian that is not a class of FT_H_csr"
    end select
end subroutine

subroutine destroy(this)
    class(FT_H_csr),intent(inout)    :: this

    if(allocated(this%csr)) deallocate(this%csr)
    if(allocated(this%i))   deallocate(this%i  )
    if(allocated(this%j))   deallocate(this%j  )
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  FEAST ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine evec_feast(this,eval,evec)
    class(FT_H_feast_csr),intent(in)       ::  this
    real(8),intent(out),allocatable        ::  eval(:)
    complex(8),intent(out),allocatable     ::  evec(:,:)

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
    fpm(3)=-nint(log10(this%ham_io%diag_acc))
    emin=this%ham_io%Ebnd(1)
    emax=this%ham_io%Ebnd(2)
    m0=this%ham_io%estNe
    if((m0.eq.0).or.(m0.gt.this%ham_io%dimH)) m0=this%ham_io%dimH
    allocate(e(m0),source=0.0d0)
    allocate(x(this%ham_io%dimh,m0),source=(0.0d0,0.0d0))
    allocate(res(this%ham_io%dimH),source=0.0d0)
    call zfeast_hcsrev('F', this%ham_io%dimH, this%csr, this%i, this%j, fpm, epsout, loop, emin, emax, m0, e, x, m, res, info)
    if(info/=0) STOP "sparse diagonalization failed"
    allocate(eval,source=e(1:m))
    deallocate(e)
    allocate(evec,source=x(1:this%ham_io%dimH,1:m))
    deallocate(x,res)
end subroutine
#endif
end module
