module m_FT_Ham_dense
!Hamiltonian type specifications using dense matrices and no external library
use m_derived_types, only: lattice, number_different_order_parameters
use m_work_ham_single, only:  work_mode, N_work
use m_H_type
use m_H_type
use m_FT_Ham_base
use m_FT_Ham_coo_rtok_base
use m_parameters_FT_Ham
use m_eigen_val_vec
use m_constants, only : pi
implicit none

type,extends(FT_Ham_base) :: FT_H_dense
    complex(8),allocatable :: Hk(:,:)             ! complex dense Hamiltonian
contains
    !necessary routines
    procedure :: init
    procedure :: set_k
    procedure :: set_work
    procedure :: calc_eval   !subroutine to calculate the eigenvalues using the internally set matrix
    procedure :: calc_evec   !subroutine to calculate the eigenvectors and eigenvalues using the internally set matrix

    !MPI
    procedure :: send
    procedure :: recv
    procedure :: distribute
    procedure :: bcast
end type


private
public FT_H_dense
contains


subroutine init(this,Hk_inp,io)
class(FT_H_dense),intent(inout)       :: this
type(H_inp_real_to_k),intent(in)      :: Hk_inp(:)
type(parameters_FT_HAM_IO),intent(in) :: io

integer :: dim_H,test

    this%dimH=3

    if (allocated(this%Hk)) ERROR STOP "H is already allocated in init FT_Ham_dense"

    allocate(this%Hk(this%dimH,this%dimH),STAT=test)
    if (test.ne.0) ERROR STOP "can not allocate Hk in init in FT_Ham_dense"

    this%Hk=(0.0d0,0.0d0)
    this%io_H=io
end subroutine

subroutine set_k(this,Hk_inp,k)
class(FT_H_dense),intent(inout)       :: this
type(H_inp_real_to_k),intent(in)      :: Hk_inp(:)
real(8),intent(in)                    :: k(3)

real(8)     :: phase_r
complex(8)  :: phase_c
integer  :: iH,i_shell,n_Ham

    this%Hk=(0.0d0,0.0d0)
    n_Ham=size(Hk_inp)
    do iH=1,n_Ham
       do i_shell=1,size(Hk_inp(iH)%H,3)
          phase_r=dot_product(Hk_inp(iH)%diffR(:,i_shell),k)
          phase_c=cmplx(cos(phase_r),sin(phase_r),8)
          this%Hk=this%Hk+phase_c*Hk_inp(iH)%H(:,:,i_shell)
       enddo
    enddo

end subroutine

subroutine set_work(this,eval,evec)
class(FT_H_dense),intent(inout)       :: this
complex(8),allocatable,intent(out)    :: eval(:)      ! array containing the eigenvalues
complex(8),allocatable,intent(out)    :: evec(:,:)   ! array containing the eigenvectors

integer :: dim_H,test

    dim_H=this%get_dimH()

    if (allocated(eval)) ERROR STOP "eigenvalues is already allocated in init FT_Ham_dense"
    if (allocated(evec)) ERROR STOP "eigenvectors is already allocated in init FT_Ham_dense"

    allocate(eval(dim_H),STAT=test,source=(0.0d0,0.0d0))
    if (test.ne.0) ERROR STOP "can not allocate eigenval in set_work in FT_Ham_dense"
    allocate(evec(dim_H,dim_H),STAT=test,source=(0.0d0,0.0d0))
    if (test.ne.0) ERROR STOP "can not allocate eigenvec in set_work in FT_Ham_dense"

end subroutine

subroutine calc_eval(this,Nin,eval,Nout)
class(FT_H_dense),intent(inout)       :: this
integer,intent(in)                    :: Nin  !size of eigenvalue input array
complex(8),intent(out)                :: eval(Nin)    !eigenvalue array
integer,intent(out)                   :: Nout !calculated number of eigenvalues

complex(8)  :: eigenvec(Nin,Nin)
real(8)     :: EPS=10d-8

eigenvec=(0.0d0,0.0d0)
Nout=3
call Jacobi(EPS,Nin,this%Hk,Nin,eval,eigenvec,Nin,1)

end subroutine

subroutine calc_evec(this,Nin,eval,evec,Nout)
class(FT_H_dense),intent(inout)     :: this
integer,intent(in)                  :: Nin  !size of eigenvalue input array
complex(8),intent(out)              :: eval(Nin)
complex(8),intent(out)              :: evec(this%dimH,Nin)
integer,intent(out)                 :: Nout !calculated number of eigenvalues

STOP 'not implemented'

end subroutine

!subroutine init_coo(this,rowind,colind,val,dim_mode,op_l,op_r,lat,mult_M_single)
!    class(FT_H_dense),intent(inout)     :: this
!    real(8),allocatable,intent(inout)   :: val(:)
!    integer,allocatable,intent(inout)   :: rowind(:),colind(:)
!    integer,intent(in)                  :: dim_mode(2)
!    character(len=*),intent(in)         :: op_l         !which order parameters are used at left  side of local Hamiltonian-matrix
!    character(len=*),intent(in)         :: op_r         !which order parameters are used at right side of local Hamiltonian-matrix
!    type(lattice),intent(in)            :: lat
!    integer,intent(in)                  :: mult_M_single !gives the multiple with which the energy_single calculation has to be multiplied (1 for on-site terms, 2 for eg. magnetic exchange)
!
!    integer :: i
!
!    allocate(this%H(this%dim_mode(1),this%dim_mode(2)),source=cmplx(0.0d0,0.0d0,8))
!    do i=1,size(val)
!        this%H(rowind(i),colind(i))=val(i)
!    enddo
!
!end subroutine

!subroutine finish_setup(this)
!    class(FT_H_dense),intent(inout)    :: this
!
!    call this%finish_setup_base()
!
!end subroutine

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
    Call MPI_Send(this%Hk, size(this%Hk), MPI_DOUBLE_COMPLEX, ithread, tag, com, ierr)
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
    allocate(this%Hk(this%dimH,this%dimH))
    Call MPI_RECV(this%Hk, size(this%Hk), MPI_DOUBLE_COMPLEX, ithread, tag, com, stat, ierr)
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
        N=shape(this%Hk)
    endif
    Call MPI_Bcast(N, 2, MPI_INTEGER, comm%mas, comm%com,ierr)
    if(.not.comm%ismas)then
        allocate(this%Hk(N(1),N(2)))
    endif
    Call MPI_Bcast(this%Hk,size(this%Hk), MPI_COMPLEX8, comm%mas, comm%com,ierr)
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

end module
