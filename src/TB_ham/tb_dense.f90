module m_H_tb_dense
use m_H_tb_base
use m_H_tb_coo
private
public H_tb_dense
type,extends(H_TB_coo_based),abstract :: H_tb_dense
    private
    complex(8),allocatable :: H(:,:)
    contains
    procedure   :: set_from_Hcoo
    procedure   :: add_child
    procedure   :: copy_child
    procedure   :: destroy_child
    procedure   :: mv
    procedure   :: mult_r
    procedure   :: get_H
    procedure   :: set_H
end type

#ifdef CPP_LAPACK
public H_zheev, H_zheevd, H_zheevr, H_zheevx
type,extends(H_tb_dense)  ::  H_zheev
    contains
    procedure   :: get_eval => eval_zheev
    procedure   :: get_evec => evec_zheev
end type

type,extends(H_tb_dense)  ::  H_zheevd
    contains
    procedure   :: get_eval => eval_zheevd
    procedure   :: get_evec => evec_zheevd
end type

type,extends(H_tb_dense)  ::  H_zheevr
    contains
    procedure   :: get_eval => eval_zheevr
    procedure   :: get_evec => evec_zheevr
end type

type,extends(H_tb_dense)  ::  H_zheevx
    contains
    procedure   :: get_eval => eval_zheevx
    procedure   :: get_evec => evec_zheevx
end type
#endif

#ifdef CPP_MKL
public H_feast_den
type,extends(H_tb_dense)  ::  H_feast_den
    contains
    procedure   :: get_evec => evec_feast
end type
#endif


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  GENERAL ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_from_Hcoo(this,H_coo)
    class(H_tb_dense),intent(inout)   :: this
    type(H_tb_coo),intent(inout)      :: H_coo
    integer                 :: nnz,i 
    complex(8),allocatable  :: val(:)
    integer,allocatable     :: rowind(:),colind(:)

    Call this%init_otherH(H_coo)
    Call H_coo%pop_par(nnz,val,rowind,colind)
    allocate(this%H(this%dimH,this%dimH),source=(0.0d0,0.0d0))
    do i=1,nnz
        this%H(rowind(i),colind(i))=this%H(rowind(i),colind(i))+val(i)
    enddo
end subroutine 

subroutine get_H(this,H)
    class(H_tb_dense),intent(inout) :: this
    complex(8),intent(inout)        :: H(this%dimH,this%dimH)

    H=this%H
end subroutine

subroutine set_H(this,H)
    class(H_tb_dense),intent(inout) :: this
    complex(8),intent(in)           :: H(this%dimH,this%dimH)

    this%H=H
end subroutine


subroutine mv(this,Hout)
    class(H_tb_dense),intent(inout) :: this
    class(H_TB),intent(inout)       :: Hout
    
    select type(Hout)
    class is(H_TB_dense)
        Call Hout%init_otherH(this)
        Call move_alloc(this%H,Hout%H)
    class default
        STOP "Cannot move H_TB_dense type to Hamiltonian that is not a class of H_TB_dense"
    end select
    Call this%destroy()
end subroutine


recursive subroutine add_child(this,H_in)
    class(H_tb_dense),intent(inout)     :: this
    class(H_tb),intent(in)              :: H_in

    class(H_tb_dense),allocatable   ::  Htmp
    type(H_tb_coo)                  ::  Hcoo_tmp

    select type(H_in)
    class is(H_tb_dense)
        this%H=this%H+H_in%H
    class is(H_TB_coo)
        Call H_in%copy(Hcoo_tmp)    !make local copy since add should not destroy the added array
        allocate(Htmp,mold=this)    
        Call Htmp%set_from_Hcoo(Hcoo_tmp)   
        Call this%add(Htmp)
        Call Htmp%destroy()
        Call Hcoo_tmp%destroy()
    class default
        ERROR STOP "Cannot add H_tb_dense type with Hamiltonian that is not a class of H_tb_dense"
    end select
end subroutine

subroutine copy_child(this,Hout)
    class(H_tb_dense),intent(in)   :: this
    class(H_TB),intent(inout)      :: Hout
    
    select type(Hout)
    class is(H_tb_dense)
        allocate(Hout%H,source=this%H)
    class default
        STOP "Cannot copy H_tb_dense type with Hamiltonian that is not a class of H_tb_dense"
    end select
end subroutine

subroutine destroy_child(this)
    class(H_tb_dense),intent(inout)    :: this

    if(allocated(this%H)) deallocate(this%H)
end subroutine

subroutine mult_r(this,vec,res,alpha,beta)
    class(H_tb_dense),intent(in)    :: this
    complex(8),intent(in   )        :: vec(this%dimH)
    complex(8),intent(inout)        :: res(this%dimH)
    complex(8),intent(in),optional     :: alpha
    complex(8),intent(in),optional     :: beta
    complex(8) :: alp, bet

    if(present(alpha))then
        alp=alpha
    else
        alp=cmplx(1.0d0,0.0d0,8)
    endif
    if(present(beta))then
        bet=beta
    else
        bet=cmplx(0.0d0,0.0d0,8)
    endif
    Call zgemv( "N", this%dimH, this%dimH, alp, this%H, this%dimH, vec, 1, bet, res, 1) 
end subroutine


#ifdef CPP_LAPACK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  ZHEEV ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eval_zheev(this,eval)
    class(H_zheev),intent(in)       ::  this
    real(8),intent(out),allocatable ::  eval(:)
    !internal
    complex(8)              :: H(this%dimH,this%dimH)
    complex(8),allocatable  :: work(:)
    complex(8)              :: tmp(1)
    integer                 :: info,lwork
    real(8)                 :: rwork(max(1,this%dimH*3-2))
    external zheev

    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVALUE AS HAMILTONIAN IS NOT SET"
    allocate(eval(this%dimH),source=0.0d0)
    H=this%H
    Call zheev('N', 'U', this%dimH, H, this%dimH, eval, tmp, -1, rwork, info) 
    if(info/=0) ERROR STOP "LAPACK ERROR"
    lwork=int(tmp(1))
    allocate(work(lwork))
    Call zheev('N', 'U', this%dimH, H, this%dimH, eval, work, lwork, rwork, info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
end subroutine

subroutine evec_zheev(this,eval,evec)
    class(H_zheev),intent(in)           ::  this
    real(8),intent(out),allocatable     ::  eval(:)
    complex(8),intent(out),allocatable  ::  evec(:,:)
    !internal
    complex(8),allocatable  :: work(:)
    complex(8)              :: tmp(1)
    integer                 :: info,lwork
    real(8)                 :: rwork(max(1,this%dimH*3-2))
    external zheev
    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVALUE AS HAMILTONIAN IS NOT SET"
    allocate(eval(this%dimH),source=0.0d0)
    allocate(evec,source=this%H)
    Call zheev('V', 'U', this%dimH, evec, this%dimH, eval, tmp, -1, rwork, info) 
    if(info/=0) ERROR STOP "LAPACK ERROR"
    lwork=int(tmp(1))
    allocate(work(lwork))
    Call zheev('V', 'U', this%dimH, evec, this%dimH, eval, work, lwork, rwork, info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  ZHEEVD ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine eval_zheevd(this,eval)
    class(H_zheevd),intent(in)      ::  this
    real(8),intent(out),allocatable ::  eval(:)
    !internal
    complex(8)              :: H(this%dimH,this%dimH)
    complex(8),allocatable  :: work(:)
    real(8),allocatable     :: rwork(:)
    integer,allocatable     :: iwork(:)
    integer ::  lwork,lrwork,liwork
    integer :: size_opt(3)
    integer                 :: info
    external zheevd
    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVECTOR AS HAMILTONIAN IS NOT SET"
    allocate(eval(this%dimH),source=0.0d0)
    H=this%H
    allocate(work(1),rwork(1),iwork(1))
    size_opt=-1
    Call zheevd('N', 'U', this%dimH, H, this%dimH, eval, work, size_opt(1), rwork, size_opt(2), iwork, size_opt(3), info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    size_opt=[int(work(1)),int(rwork(1)),iwork]
    deallocate(work,rwork,iwork)
    allocate(work(size_opt(1)),rwork(size_opt(2)),iwork(size_opt(3)))
    !destroys lower triangle of H
    Call zheevd('N', 'U', this%dimH, H, this%dimH, eval, work, size_opt(1), rwork, size_opt(2), iwork, size_opt(3), info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
end subroutine

subroutine evec_zheevd(this,eval,evec)
    class(H_zheevd),intent(in)          ::  this
    real(8),intent(out),allocatable     ::  eval(:)
    complex(8),intent(out),allocatable  ::  evec(:,:)
    !internal
    complex(8),allocatable  :: work(:)
    real(8),allocatable     :: rwork(:)
    integer,allocatable     :: iwork(:)
    integer :: size_opt(3)
    integer                 :: info
    external zheevd
    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVECTOR AS HAMILTONIAN IS NOT SET"
    allocate(eval(this%dimH),source=0.0d0)
    allocate(evec,source=this%H)
    allocate(work(1),rwork(1),iwork(1))
    size_opt=-1
    Call zheevd('V', 'U', this%dimH, evec, this%dimH, eval, work, size_opt(1), rwork, size_opt(2), iwork, size_opt(3), info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    size_opt=[int(work(1)),int(rwork(1)),iwork]
    deallocate(work,rwork,iwork)
    allocate(work(size_opt(1)),rwork(size_opt(2)),iwork(size_opt(3)))
    !destroys lower triangle of H
    Call zheevd('V', 'U', this%dimH, evec, this%dimH, eval, work, size_opt(1), rwork, size_opt(2), iwork, size_opt(3), info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  ZHEEVR ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine eval_zheevr(this,eval)
    class(H_zheevr),intent(in)      ::  this
    real(8),intent(out),allocatable ::  eval(:)
    !internal
    complex(8)              :: H(this%dimH,this%dimH)
    real(8)                 :: w(this%dimH) !eigenvalues
    !temporary Lapack values
    complex(8),allocatable  :: work(:)
    real(8),allocatable     :: rwork(:)
    integer,allocatable     :: iwork(:)
    integer :: size_opt(3)
    integer                 :: isuppz(2*this%dimH)
    real(8)                 :: abstol
    integer                 :: info,lwork
    integer                 :: Nev  !number eigenvalues found
    complex(8)              :: z(1) !not referenced here
    external zheevr
    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVALUE AS HAMILTONIAN IS NOT SET"
    H=this%H
    abstol=this%diag_acc
    allocate(work(1),rwork(1),iwork(1))
    size_opt=-1
    Call zheevr('N', 'V', 'U', this%dimH, H, this%dimH, this%Ebnd(1), this%Ebnd(2), 0, 0, abstol, Nev, w, z, 1, isuppz, & 
                & work, size_opt(1), rwork, size_opt(2), iwork, size_opt(3), info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    size_opt=[int(work(1)),int(rwork(1)),iwork]
    deallocate(work,rwork,iwork)
    allocate(work(size_opt(1)),rwork(size_opt(2)),iwork(size_opt(3)))
    Call zheevr('N', 'V', 'U', this%dimH, H, this%dimH, this%Ebnd(1), this%Ebnd(2), 0, 0, abstol, Nev, w, z, 1, isuppz, & 
                & work, size_opt(1), rwork, size_opt(2), iwork, size_opt(3), info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    allocate(eval,source=w(1:Nev))
end subroutine


subroutine evec_zheevr(this,eval,evec)
    class(H_zheevr),intent(in)          ::  this
    real(8),intent(out),allocatable     ::  eval(:)
    complex(8),intent(out),allocatable  ::  evec(:,:)
    !internal
    complex(8)              :: H(this%dimH,this%dimH)
    real(8)                 :: w(this%dimH) !eigenvalues
    !temporary Lapack values
    complex(8),allocatable  :: work(:)
    real(8),allocatable     :: rwork(:)
    integer,allocatable     :: iwork(:)
    integer :: size_opt(3)
    integer                 :: isuppz(2*this%dimH)
    real(8)                 :: abstol
    integer                 :: info,lwork
    integer                 :: Nev  !number eigenvalues found
    complex(8),allocatable  :: z(:,:) !internal eigenvectors
    external zheevr

    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVALUE AS HAMILTONIAN IS NOT SET"
    H=this%H
    abstol=this%diag_acc
    allocate(work(1),rwork(1),iwork(1))
    allocate(z(this%dimH,this%dimH))   !reduces memory usage with this?
    size_opt=-1
    Call zheevr('V', 'V', 'U', this%dimH, H, this%dimH, this%Ebnd(1), this%Ebnd(2), 0, 0, abstol, Nev, w, z, this%dimH, isuppz, & 
                & work, size_opt(1), rwork, size_opt(2), iwork, size_opt(3), info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    size_opt=[int(work(1)),int(rwork(1)),iwork]
    deallocate(work,rwork,iwork)
    allocate(work(size_opt(1)),rwork(size_opt(2)),iwork(size_opt(3)))
    Call zheevr('V', 'V', 'U', this%dimH, H, this%dimH, this%Ebnd(1), this%Ebnd(2), 0, 0, abstol, Nev, w, z, this%dimH, isuppz, & 
                & work, size_opt(1), rwork, size_opt(2), iwork, size_opt(3), info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    allocate(eval,source=w(1:Nev))
    allocate(evec,source=z(1:this%dimH,1:Nev))
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  ZHEEVX ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine eval_zheevx(this,eval)
    class(H_zheevx),intent(in)      ::  this
    real(8),intent(out),allocatable ::  eval(:)
    !internal
    complex(8)              :: H(this%dimH,this%dimH)
    real(8)                 :: w(this%dimH) !eigenvalues
    !temporary Lapack values
    complex(8),allocatable  :: work(:)
    real(8)                 :: rwork(7*this%dimH)
    integer                 :: iwork(5*this%dimH)
    integer                 :: ifail(this%dimH)
    real(8)                 :: abstol
    integer                 :: info,lwork
    integer                 :: Nev  !number eigenvalues found
    complex(8)              :: z(1) !not referenced here
    complex(8)              :: work_tmp(1)
    external zheevx
    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVALUE AS HAMILTONIAN IS NOT SET"
    H=this%H
    abstol=this%diag_acc
    Call zheevx('N', 'V', 'U', this%dimH, H, this%dimH, this%Ebnd(1), this%Ebnd(2), 0, 0, abstol, Nev, w, z, 1, & 
                & work_tmp, -1, rwork, iwork, ifail, info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    lwork=int(work_tmp(1))
    allocate(work(lwork))
    Call zheevx('N', 'V', 'U', this%dimH, H, this%dimH, this%Ebnd(1), this%Ebnd(2), 0, 0, abstol, Nev, w, z, 1, & 
                & work, lwork, rwork, iwork, ifail, info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    allocate(eval,source=w(1:Nev))
end subroutine


subroutine evec_zheevx(this,eval,evec)
    class(H_zheevx),intent(in)          ::  this
    real(8),intent(out),allocatable     ::  eval(:)
    complex(8),intent(out),allocatable  ::  evec(:,:)
    !internal
    complex(8)              :: H(this%dimH,this%dimH)
    real(8)                 :: w(this%dimH) !eigenvalues
    complex(8),allocatable  :: work(:)
    real(8)                 :: rwork(7*this%dimH)
    integer                 :: iwork(5*this%dimH)
    integer                 :: ifail(this%dimH)
    real(8)                 :: abstol
    integer                 :: info,lwork
    integer                 :: Nev  !number eigenvalues found
    complex(8),allocatable  :: z(:,:)
    complex(8)              :: work_tmp(1)
    external zheevx

    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVALUE AS HAMILTONIAN IS NOT SET"
    H=this%H
    abstol=this%diag_acc
    allocate(z(this%dimH,this%dimH))   !reduces memory usage with this?
    Call zheevx('V', 'V', 'U', this%dimH, H, this%dimH, this%Ebnd(1), this%Ebnd(2), 0, 0, abstol, Nev, w, z, this%dimH,& 
                & work_tmp, -1, rwork, iwork, ifail, info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    lwork=int(work_tmp(1))
    allocate(work(lwork))
    Call zheevx('V', 'V', 'U', this%dimH, H, this%dimH, this%Ebnd(1), this%Ebnd(2), 0, 0, abstol, Nev, w, z, this%dimH, & 
                & work, lwork, rwork, iwork, ifail, info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    allocate(eval,source=w(1:Nev))
    allocate(evec,source=z(1:this%dimH,1:Nev))
end subroutine


#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  FEAST ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef CPP_MKL
subroutine evec_feast(this,eval,evec)
    class(H_feast_den),intent(in)       ::  this
    real(8),intent(out),allocatable     ::  eval(:)
    complex(8),intent(out),allocatable  ::  evec(:,:)

    integer                     :: fpm(128)
    real(8)                     :: epsout
    integer                     :: loop
    integer                     :: m0,m
    real(8),allocatable         :: res(:)
    integer                     :: info
    complex(8),allocatable      :: x(:,:)
    real(8),allocatable         :: e(:)

    Call feastinit(fpm) 
    fpm(1)=1
    fpm(2)=8
    fpm(3)=-nint(log10(this%diag_acc))
    m0=this%estNe
    if(m0==0.or.m0>this%dimH) m0=this%dimH
    allocate(res(m0),e(m0),x(this%dimH,m0))
    call zfeast_heev ( 'U', this%dimH, this%H, this%dimH, fpm, epsout, loop, this%Ebnd(1), this%Ebnd(2), m0, e, x, m, res, info )
    if(info/=0) STOP 'info of zfest_heev not zero'
    allocate(eval,source=e(1:m))
    allocate(evec,source=x(1:this%dimH,1:m))
    deallocate(x,e)
end subroutine
#endif

end module
