module m_FT_Ham_dense
!Hamiltonian type specifications using dense matrices and no external library
use m_derived_types, only: lattice, number_different_order_parameters
use m_work_ham_single, only:  work_mode, N_work
use m_H_type
use m_parameters_FT_Ham
use m_H_type
use m_FT_Ham_coo
implicit none

type,extends(t_H_base),abstract :: FT_H_dense
    complex(8),allocatable       :: H(:,:)
    type(parameters_FT_HAM_IO)   :: io_H
contains
    !necessary t_H routines
    procedure :: set_from_Hcoo
    procedure :: init_coo
    procedure :: init_connect
    procedure :: init_mult_connect_2
    procedure :: finish_setup

    procedure :: add_child
    procedure :: destroy_child
    procedure :: copy_child

    procedure :: optimize
    procedure :: mult_r,mult_l
    !MPI
    procedure :: send
    procedure :: recv
    procedure :: distribute
    procedure :: bcast
end type

#ifdef CPP_LAPACK
public H_zheev, H_zheevd, H_zheevr
type,extends(FT_H_dense)  ::  H_zheev
    contains
    procedure   :: get_eval => eval_zheev
    procedure   :: get_evec => evec_zheev
end type

type,extends(FT_H_dense)  ::  H_zheevd
    contains
    procedure   :: get_eval => eval_zheevd
    procedure   :: get_evec => evec_zheevd
end type

type,extends(FT_H_dense)  ::  H_zheevr
    contains
    procedure   :: get_eval => eval_zheevr
    procedure   :: get_evec => evec_zheevr
end type
#endif

!#ifdef CPP_MKL
!public H_feast_den
!type,extends(FT_H_dense)  ::  H_feast_den
!    contains
!    procedure   :: get_evec => evec_feast
!end type
!#endif

private
public FT_H_dense
contains

subroutine mult_r(this,lat,res,work,alpha,beta)
    !mult
    use m_derived_types, only: lattice
    class(FT_h_dense),intent(in)     :: this
    type(lattice), intent(in)        :: lat
    real(8), intent(inout)           :: res(:)   !result matrix-vector product
    type(work_mode),intent(inout)    :: work
    real(8),intent(in),optional      :: alpha
    real(8),intent(in),optional      :: beta
    ! internal

    STOP "Cannot use the mult_r with FT_h_dense type with Hamiltonian because type is CMPLX"

end subroutine

subroutine mult_l(this,lat,res,work,alpha,beta)
    use m_derived_types, only: lattice
    class(FT_H_dense),intent(in)     :: this
    type(lattice), intent(in)        :: lat
    real(8), intent(inout)           :: res(:)
    type(work_mode),intent(inout)    :: work
    real(8),intent(in),optional      :: alpha
    real(8),intent(in),optional      :: beta
    ! internal

    STOP "Cannot use the mult_l with FT_h_dense type with Hamiltonian because type is CMPLX"

end subroutine


subroutine optimize(this)
    class(FT_H_dense),intent(inout)   :: this

    !nothing to optimize here
    continue
end subroutine


subroutine copy_child(this,Hout)
    class(FT_H_dense),intent(in)         :: this
    class(t_H_base),intent(inout)        :: Hout

    select type(Hout)
    class is(H_inp_k_coo)
        Call this%copy_child(Hout)
    class default
        STOP "Cannot copy FT_h_dense type with Hamiltonian that is not a class of t_H_base"
    end select
end subroutine


subroutine add_child(this,H_in)
    class(FT_H_dense),intent(inout)    :: this
    class(t_H_base),intent(in)         :: H_in

    integer :: i

    select type(H_in)
    class is(H_inp_k_coo)
        call this%add_child(H_in)
    class default
        STOP "Cannot add FT_h_dense type with Hamiltonian that is not a class of t_H_base"
    end select

end subroutine

subroutine destroy_child(this)
    class(FT_H_dense),intent(inout)    :: this

    if(this%is_set())then
        deallocate(this%H)
    endif
end subroutine

subroutine set_from_Hcoo(this,H_coo)
    class(FT_H_dense),intent(inout)  :: this
    type(H_inp_k_coo),intent(inout)  :: H_coo

    !local
    integer                 :: nnz,i
    complex(8),allocatable  :: val(:)
    integer,allocatable     :: rowind(:),colind(:)

    Call H_coo%pop_par(this%dimH,nnz,val,rowind,colind)
    allocate(this%H(this%dimH(1),this%dimH(2)),source=cmplx(0.0d0,0.0d0,8))
    do i=1,nnz
        this%H(rowind(i),colind(i))=val(i)
    enddo
end subroutine

subroutine init_coo(this,rowind,colind,val,dim_mode,op_l,op_r,lat,mult_M_single)
    class(FT_H_dense),intent(inout)     :: this
    real(8),allocatable,intent(inout)   :: val(:)
    integer,allocatable,intent(inout)   :: rowind(:),colind(:)
    integer,intent(in)                  :: dim_mode(2)
    character(len=*),intent(in)         :: op_l         !which order parameters are used at left  side of local Hamiltonian-matrix
    character(len=*),intent(in)         :: op_r         !which order parameters are used at right side of local Hamiltonian-matrix
    type(lattice),intent(in)            :: lat
    integer,intent(in)                  :: mult_M_single !gives the multiple with which the energy_single calculation has to be multiplied (1 for on-site terms, 2 for eg. magnetic exchange)

    integer :: i

    allocate(this%H(this%dim_mode(1),this%dim_mode(2)),source=cmplx(0.0d0,0.0d0,8))
    do i=1,size(val)
        this%H(rowind(i),colind(i))=val(i)
    enddo

end subroutine

subroutine init_connect(this,connect,Hval,Hval_ind,order,lat,mult_M_single)
    class(FT_H_dense),intent(inout)   :: this
    type(lattice),intent(in)        :: lat
    character(2),intent(in)         :: order
    real(8),intent(in)              :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
    integer,intent(in)              :: Hval_ind(:,:)
    integer,intent(in)              :: connect(:,:)
    integer,intent(in)              :: mult_M_single

    STOP "Cannot use init_connect with FT_h_dense type with Hamiltonian because type is CMPLX"

end subroutine

subroutine init_mult_connect_2(this,connect,Hval,Hval_ind,op_l,op_r,lat,mult_M_single,dim_mode_in)
    class(FT_H_dense),intent(inout)   :: this
    type(lattice),intent(in)        :: lat
    character(*),intent(in)         :: op_l
    character(*),intent(in)         :: op_r
    real(8),intent(in)              :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
    integer,intent(in)              :: Hval_ind(:,:)
    integer,intent(in)              :: connect(:,:)
    integer,intent(in)              :: mult_M_single
    integer,intent(in),optional     :: dim_mode_in(2)   !optional way of putting in dim_mode directly (mainly for custom(not fully unfolded)rankN tensors)

    STOP "Cannot use init_mult_connect_2 with FT_h_dense type with Hamiltonian because type is CMPLX"

end subroutine

subroutine finish_setup(this)
    class(FT_H_dense),intent(inout)    :: this

    call this%finish_setup_base()

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!            MPI ROUTINES           !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine send(this,ithread,tag,com)
    use mpi_basic
    class(FT_h_dense),intent(in)     :: this
    integer,intent(in)               :: ithread
    integer,intent(in)               :: tag
    integer,intent(in)               :: com

#ifdef CPP_MPI
    integer     :: ierr
    Call this%send_base(ithread,tag,com)
    Call MPI_Send(this%H, size(this%H), MPI_DOUBLE_COMPLEX, ithread, tag, com, ierr)
#else
    continue
#endif
end subroutine

subroutine recv(this,ithread,tag,com)
    use mpi_basic
    class(FT_h_dense),intent(inout)  :: this
    integer,intent(in)               :: ithread
    integer,intent(in)               :: tag
    integer,intent(in)               :: com

#ifdef CPP_MPI
    integer     :: stat(MPI_STATUS_SIZE)
    integer     :: ierr

    Call this%recv_base(ithread,tag,com)
    allocate(this%H(this%dimH(1),this%dimH(2)))
    Call MPI_RECV(this%H, size(this%H), MPI_DOUBLE_COMPLEX, ithread, tag, com, stat, ierr)
#else
    continue
#endif
end subroutine

subroutine bcast(this,comm)
    use mpi_basic
    class(FT_h_dense),intent(inout)   ::  this
    type(mpi_type),intent(in)         ::  comm
#ifdef CPP_MPI
    integer     :: ierr
    integer     :: N(2)

    Call this%bcast_base(comm)
    if(comm%ismas)then
        N=shape(this%H)
    endif
    Call MPI_Bcast(N, 2, MPI_INTEGER, comm%mas, comm%com,ierr)
    if(.not.comm%ismas)then
        allocate(this%H(N(1),N(2)))
    endif
    Call MPI_Bcast(this%H,size(this%H), MPI_COMPLEX8, comm%mas, comm%com,ierr)
#else
    continue
#endif
end subroutine

subroutine distribute(this,comm)
    use mpi_basic
    class(FT_h_dense),intent(inout)  ::  this
    type(mpi_type),intent(in)        ::  comm

    ERROR STOP "Inner MPI-parallelization of the FT Hamiltonian is not implemented for dense Hamiltonians"
end subroutine

#ifdef CPP_LAPACK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  ZHEEV ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eval_zheev(this,eval)
    class(H_zheev),intent(in)       ::  this
    real(8),intent(out),allocatable ::  eval(:)
    !internal
    complex(8)              :: H(this%io_H%dimH,this%io_H%dimH)
    complex(8),allocatable  :: work(:)
    complex(8)              :: tmp(1)
    integer                 :: info,lwork
    real(8)                 :: rwork(max(1,this%io_H%dimH*3-2))
    external zheev

    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVALUE AS HAMILTONIAN IS NOT SET"
    allocate(eval(this%io_H%dimH),source=0.0d0)
    H=this%H
    Call zheev('N', 'U', this%io_H%dimH, H, this%io_H%dimH, eval, tmp, -1, rwork, info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    lwork=int(tmp(1))
    allocate(work(lwork))
    Call zheev('N', 'U', this%io_H%dimH, H, this%io_H%dimH, eval, work, lwork, rwork, info)
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
    real(8)                 :: rwork(max(1,this%io_H%dimH*3-2))
    external zheev
    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVALUE AS HAMILTONIAN IS NOT SET"
    allocate(eval(this%io_H%dimH),source=0.0d0)
    allocate(evec,source=this%H)
    Call zheev('V', 'U', this%io_H%dimH, evec, this%io_H%dimH, eval, tmp, -1, rwork, info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    lwork=int(tmp(1))
    allocate(work(lwork))
    Call zheev('V', 'U', this%io_H%dimH, evec, this%io_H%dimH, eval, work, lwork, rwork, info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  ZHEEVD ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine eval_zheevd(this,eval)
    class(H_zheevd),intent(in)      ::  this
    real(8),intent(out),allocatable ::  eval(:)
    !internal
    complex(8)              :: H(this%io_H%dimH,this%io_H%dimH)
    complex(8),allocatable  :: work(:)
    real(8),allocatable     :: rwork(:)
    integer,allocatable     :: iwork(:)
    integer ::  lwork,lrwork,liwork
    integer :: size_opt(3)
    integer                 :: info
    external zheevd
    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVECTOR AS HAMILTONIAN IS NOT SET"
    allocate(eval(this%io_H%dimH),source=0.0d0)
    H=this%H
    allocate(work(1),rwork(1),iwork(1))
    size_opt=-1
    Call zheevd('N', 'U', this%io_H%dimH, H, this%io_H%dimH, eval, work, size_opt(1), rwork, size_opt(2), iwork, size_opt(3), info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    size_opt=[int(work(1)),int(rwork(1)),iwork]
    deallocate(work,rwork,iwork)
    allocate(work(size_opt(1)),rwork(size_opt(2)),iwork(size_opt(3)))
    !destroys lower triangle of H
    Call zheevd('N', 'U', this%io_H%dimH, H, this%io_H%dimH, eval, work, size_opt(1), rwork, size_opt(2), iwork, size_opt(3), info)
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
    allocate(eval(this%io_H%dimH),source=0.0d0)
    allocate(evec,source=this%H)
    allocate(work(1),rwork(1),iwork(1))
    size_opt=-1
    Call zheevd('V', 'U', this%io_H%dimH, evec, this%io_H%dimH, eval, work, size_opt(1), rwork, size_opt(2), iwork, size_opt(3), info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    size_opt=[int(work(1)),int(rwork(1)),iwork]
    deallocate(work,rwork,iwork)
    allocate(work(size_opt(1)),rwork(size_opt(2)),iwork(size_opt(3)))
    !destroys lower triangle of H
    Call zheevd('V', 'U', this%io_H%dimH, evec, this%io_H%dimH, eval, work, size_opt(1), rwork, size_opt(2), iwork, size_opt(3), info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  ZHEEVR ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine eval_zheevr(this,eval)
    class(H_zheevr),intent(in)      ::  this
    real(8),intent(out),allocatable ::  eval(:)
    !internal
    complex(8)              :: H(this%io_H%dimH,this%io_H%dimH)
    real(8)                 :: w(this%io_H%dimH) !eigenvalues
    !temporary Lapack values
    complex(8),allocatable  :: work(:)
    real(8),allocatable     :: rwork(:)
    integer,allocatable     :: iwork(:)
    integer :: size_opt(3)
    integer                 :: isuppz(2*this%io_H%dimH)
    real(8)                 :: abstol
    integer                 :: info,lwork
    integer                 :: Nev  !number eigenvalues found
    complex(8)              :: z(1) !not referenced here
    external zheevr
    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVALUE AS HAMILTONIAN IS NOT SET"
    H=this%H
    abstol=this%io_H%diag_acc
    allocate(work(1),rwork(1),iwork(1))
    size_opt=-1
    Call zheevr('N', 'V', 'U', this%io_H%dimH, H, this%io_H%dimH, this%io_H%Ebnd(1), this%io_H%Ebnd(2), 0, 0, abstol, Nev, w, z, 1, isuppz, &
                & work, size_opt(1), rwork, size_opt(2), iwork, size_opt(3), info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    size_opt=[int(work(1)),int(rwork(1)),iwork]
    deallocate(work,rwork,iwork)
    allocate(work(size_opt(1)),rwork(size_opt(2)),iwork(size_opt(3)))
    Call zheevr('N', 'V', 'U', this%io_H%dimH, H, this%io_H%dimH, this%io_H%Ebnd(1), this%io_H%Ebnd(2), 0, 0, abstol, Nev, w, z, 1, isuppz, &
                & work, size_opt(1), rwork, size_opt(2), iwork, size_opt(3), info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    allocate(eval,source=w(1:Nev))
end subroutine


subroutine evec_zheevr(this,eval,evec)
    class(H_zheevr),intent(in)          ::  this
    real(8),intent(out),allocatable     ::  eval(:)
    complex(8),intent(out),allocatable  ::  evec(:,:)
    !internal
    complex(8)              :: H(this%io_H%dimH,this%io_H%dimH)
    real(8)                 :: w(this%io_H%dimH) !eigenvalues
    !temporary Lapack values
    complex(8),allocatable  :: work(:)
    real(8),allocatable     :: rwork(:)
    integer,allocatable     :: iwork(:)
    integer :: size_opt(3)
    integer                 :: isuppz(2*this%io_H%dimH)
    real(8)                 :: abstol
    integer                 :: info,lwork
    integer                 :: Nev  !number eigenvalues found
    complex(8),allocatable  :: z(:,:) !internal eigenvectors
    external zheevr

    if(.not.this%is_set()) ERROR STOP "CANNOT GET EIGENVALUE AS HAMILTONIAN IS NOT SET"
    H=this%H
    abstol=this%io_H%diag_acc
    allocate(work(1),rwork(1),iwork(1))
!    if(this%estNe<1) STOP "estimated number of eigenvalues must be at least 1('TB_diag_estNe')"
!    allocate(z(this%dimH,min(this%estNe,this%dimH)))   !reduces memory usage with this?
    allocate(z(this%io_H%dimH,this%io_H%dimH))   !reduces memory usage with this?
    size_opt=-1
    Call zheevr('V', 'V', 'U', this%io_H%dimH, H, this%io_H%dimH, this%io_H%Ebnd(1), this%io_H%Ebnd(2), 0, 0, abstol, Nev, w, z, this%io_H%dimH, isuppz, &
                & work, size_opt(1), rwork, size_opt(2), iwork, size_opt(3), info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
    size_opt=[int(work(1)),int(rwork(1)),iwork]
    deallocate(work,rwork,iwork)
    allocate(work(size_opt(1)),rwork(size_opt(2)),iwork(size_opt(3)))
!    write(*,*) "If this crashes, increase TB_diag_estNe"
    Call zheevr('V', 'V', 'U', this%io_H%dimH, H, this%io_H%dimH, this%io_H%Ebnd(1), this%io_H%Ebnd(2), 0, 0, abstol, Nev, w, z, this%io_H%dimH, isuppz, &
                & work, size_opt(1), rwork, size_opt(2), iwork, size_opt(3), info)
    if(info/=0) ERROR STOP "LAPACK ERROR"
!    write(*,*) "Did not crash"
    allocate(eval,source=w(1:Nev))
    allocate(evec,source=z(1:this%io_H%dimH,1:Nev))
end subroutine

#endif
end module
