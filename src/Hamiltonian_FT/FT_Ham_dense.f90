module m_FT_Ham_dense
!Hamiltonian type specifications using dense matrices and no external library
use m_derived_types, only: lattice, number_different_order_parameters
use m_work_ham_single, only:  work_mode, N_work
use m_H_type
use m_H_type
use m_FT_Ham_base
use m_FT_Ham_coo_rtok_base
implicit none

type,extends(FT_Ham_base),abstract :: FT_H_dense
    complex(8),allocatable          :: H(:,:)         ! complex dense Hamiltonian
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
type(H_inp_real_to_k),intent(inout)   :: Hk_inp
type(parameters_TB_IO_H),intent(in)   :: io

integer :: dim_H,test
    dim_H=this%get_dimH()

    if allocated(this%H) ERROR STOP "H is already allocated in init FT_Ham_dense"

    allocate(this%H(dim_H,dim_H),STAT=test)
write(*,*) test

pause
this%H=(0.0d0,0.0d0)
end subroutine

subroutine set_k(this,Hr,k)
class(FT_H_dense),intent(inout)       :: this
real(8),intent(in)                    :: Hr(:,:,:)
real(8),intent(in)                    :: k(3)

real(8)     :: phase_r
complex(8)  :: phase_c
integer  :: iH

    do iH=1,size(Hr,3)
       phase_r=dot_product(this%diffR(:,iH),k)
       phase_c=cmplx(cos(phase_r),sin(phase_r),8)
       this%H=this%H+phase_c*this%H_R(:,:,iH)
    enddo

end subroutine

subroutine set_work(this,Hr,k)
class(FT_H_dense),intent(inout)       :: this
real(8),intent(in)                    :: Hr(:,:,:)
real(8),intent(in)                    :: k(3)

real(8)     :: phase_r
complex(8)  :: phase_c
integer  :: iH
    class(FT_H_dense),intent(inout)      :: this
    type(work_ham),intent(inout)        :: work
    integer                             :: sizes(N_work)

    sizes=0
    sizes(1)=max(1,3*this%dimH-2)   !real(RWORK)
    sizes(3)=max(1,2*this%dimH-1)   !complex(WORK)
    Call work%set(sizes)

end subroutine

!subroutine copy_child(this,Hout)
!    class(FT_H_dense),intent(in)         :: this
!    class(t_H_base),intent(inout)        :: Hout
!
!    select type(Hout)
!    class is(H_inp_k_coo)
!        Call this%copy_child(Hout)
!    class default
!        STOP "Cannot copy FT_h_dense type with Hamiltonian that is not a class of t_H_base"
!    end select
!end subroutine


!subroutine add_child(this,H_in)
!    class(FT_H_dense),intent(inout)    :: this
!    class(t_H_base),intent(in)         :: H_in
!
!    integer :: i
!
!    select type(H_in)
!    class is(H_inp_k_coo)
!        call this%add_child(H_in)
!    class default
!        STOP "Cannot add FT_h_dense type with Hamiltonian that is not a class of t_H_base"
!    end select
!
!end subroutine

!subroutine destroy_child(this)
!    class(FT_H_dense),intent(inout)    :: this
!
!    if(this%is_set())then
!        deallocate(this%H)
!    endif
!end subroutine

!subroutine set_from_Hcoo(this,H_coo)
!    class(FT_H_dense),intent(inout)  :: this
!    type(H_inp_k_coo),intent(inout)  :: H_coo
!
!    !local
!    integer                 :: nnz,i
!    complex(8),allocatable  :: val(:)
!    integer,allocatable     :: rowind(:),colind(:)
!
!    Call H_coo%pop_par(this%dimH,nnz,val,rowind,colind)
!    allocate(this%H(this%dimH(1),this%dimH(2)),source=cmplx(0.0d0,0.0d0,8))
!    do i=1,nnz
!        this%H(rowind(i),colind(i))=val(i)
!    enddo
!end subroutine

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

end module
