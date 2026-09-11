module m_H_dense
!Hamiltonian type specifications using dense matrices and no external library
use m_derived_types, only: lattice, number_different_order_parameters
use m_H_coo_based
use m_work_ham_single, only:  work_mode, N_work

type,extends(t_H_coo_based) :: t_H_dense
    real(8),allocatable   :: H(:,:)
contains
    !necessary t_H routines
    procedure :: set_from_Hcoo 

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

private
public t_H,t_H_dense
contains 

subroutine mult_r(this,lat,res,work,alpha,beta)
    !mult
    use m_derived_types, only: lattice
    class(t_h_dense),intent(in)     :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)   !result matrix-vector product
    type(work_mode),intent(inout)   :: work
    real(8),intent(in),optional     :: alpha
    real(8),intent(in),optional     :: beta
    ! internal
    real(8),pointer,contiguous      :: modes(:)
    integer                         :: work_size(N_work)

    if(size(res)/=this%dimH(1)) STOP "size of vec is wrong"

    Call this%mode_r%get_mode(lat,modes,work,work_size)
    if(present(alpha).and.present(beta))then
        res=alpha*matmul(this%H,modes)+beta*res
    elseif(present(alpha))then
        res=alpha*matmul(this%H,modes)
    elseif(present(beta))then
        res=matmul(this%H,modes)+beta*res
    else
        res=matmul(this%H,modes)
    endif

    nullify(modes)
    work%offset=work%offset-work_size

end subroutine 

subroutine mult_l(this,lat,res,work,alpha,beta)
    use m_derived_types, only: lattice
    class(t_h_dense),intent(in)     :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)
    type(work_mode),intent(inout)   :: work
    real(8),intent(in),optional     :: alpha
    real(8),intent(in),optional     :: beta
    ! internal
    real(8),pointer,contiguous      :: modes(:)
    integer                         :: work_size(N_work)

    if(size(res)/=this%dimH(2)) STOP "size of vec is wrong"

    Call this%mode_l%get_mode(lat,modes,work,work_size)
    if(present(alpha).and.present(beta))then
        res=alpha*matmul(modes,this%H)+beta*res
    elseif(present(alpha))then
        res=alpha*matmul(modes,this%H)+res
    elseif(present(beta))then
        res=matmul(modes,this%H)+beta*res
    else
        res=matmul(modes,this%H)
    endif
    nullify(modes)
    work%offset=work%offset-work_size
end subroutine 


subroutine optimize(this)
    class(t_h_dense),intent(inout)   :: this

    !nothing to optimize here
    continue 
end subroutine


subroutine copy_child(this,Hout)
    class(t_h_dense),intent(in)   :: this
    class(t_H_base),intent(inout)        :: Hout
    
    select type(Hout)
    class is(t_h_dense)
        allocate(Hout%H,source=this%H)
        Call this%copy_deriv(Hout)
    class default
        STOP "Cannot copy t_h_dense type with Hamiltonian that is not a class of t_h_dense"
    end select
end subroutine


subroutine add_child(this,H_in)
    class(t_h_dense),intent(inout)    :: this
    class(t_H_base),intent(in)             :: H_in

    select type(H_in)
    class is(t_h_dense)
        this%H=this%H+H_in%H
    class default
        STOP "Cannot add t_h_dense type with Hamiltonian that is not a class of t_h_dense"
    end select

end subroutine 

subroutine destroy_child(this)
    class(t_h_dense),intent(inout)    :: this

    if(this%is_set())then
        deallocate(this%H)
    endif
end subroutine

subroutine set_from_Hcoo(this,H_coo)
    class(t_h_dense),intent(inout)  :: this
    type(t_H_coo),intent(inout)     :: H_coo

    !local
    integer                 :: nnz,i
    real(8),allocatable     :: val(:)
    integer,allocatable     :: rowind(:),colind(:)

    Call H_coo%pop_par(this%dimH,nnz,val,rowind,colind)
    allocate(this%H(this%dimH(1),this%dimH(2)),source=0.0d0)
    do i=1,nnz
        this%H(rowind(i),colind(i))=val(i)
    enddo
end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!            MPI ROUTINES           !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine send(this,ithread,tag,com)
    use mpi_basic                
    class(t_H_dense),intent(in)     :: this
    integer,intent(in)              :: ithread
    integer,intent(in)              :: tag
    integer,intent(in)              :: com

#ifdef CPP_MPI
    integer     :: ierr
    Call this%send_base(ithread,tag,com)
    Call MPI_Send(this%H, size(this%H), MPI_DOUBLE_PRECISION, ithread, tag, com, ierr)
#else
    continue
#endif
end subroutine

subroutine recv(this,ithread,tag,com)
    use mpi_basic                
    class(t_H_dense),intent(inout)  :: this
    integer,intent(in)              :: ithread
    integer,intent(in)              :: tag
    integer,intent(in)              :: com

#ifdef CPP_MPI
    integer     :: stat(MPI_STATUS_SIZE) 
    integer     :: ierr
    
    Call this%recv_base(ithread,tag,com)
    allocate(this%H(this%dimH(1),this%dimH(2)))
    Call MPI_RECV(this%H, size(this%H), MPI_DOUBLE_PRECISION, ithread, tag, com, stat, ierr)
#else
    continue
#endif
end subroutine

subroutine bcast(this,comm)
    use mpi_basic                
    class(t_H_dense),intent(inout)        ::  this
    type(mpi_type),intent(in)       ::  comm
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
    Call MPI_Bcast(this%H,size(this%H), MPI_REAL8, comm%mas, comm%com,ierr)
#else
    continue
#endif
end subroutine 

subroutine distribute(this,comm)
    use mpi_basic                
    class(t_H_dense),intent(inout)  ::  this
    type(mpi_type),intent(in)       ::  comm

    ERROR STOP "Inner MPI-parallelization of the Hamiltonian is not implemented for dense Hamiltonians"
end subroutine 

end module
