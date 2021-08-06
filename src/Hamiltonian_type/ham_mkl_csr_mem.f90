module m_H_mkl_csr_mem
#if defined(CPP_MKL)
!Hamiltonian type specifications using MKL_SPARSE inspector mkl in csr 
use MKL_SPBLAS
use m_H_mkl_csr
use m_type_lattice, only: dim_modes_inner, lattice,number_different_order_parameters
use m_H_coo_based
use mkl_spblas_util, only: unpack_csr
use m_work_ham_single
USE, INTRINSIC :: ISO_C_BINDING 

private
public t_H,t_H_mkl_csr_mem

type,extends(t_H_mkl_csr) :: t_H_mkl_csr_mem
    type(SPARSE_MATRIX_T)   :: HT
    !type(matrix_descr)      :: descr !should be same descr as in this%H

    !pointers to Hamiltonian data handled by mkl (row major format)
    real(C_DOUBLE),pointer  :: HT_val(:)
    integer(C_INT),pointer  :: HT_inner(:),HT_outer(:)

contains
    procedure :: set_from_Hcoo

    procedure :: add_child 
    procedure :: destroy_child    
    procedure :: copy_child 

    procedure :: optimize

    procedure :: mult_l
    procedure :: mult_l_disc

    procedure :: set_work_single
    procedure :: get_work_size_single

    procedure :: set_auxiliaries

    !MPI
    procedure :: send
    procedure :: recv
    procedure :: distribute
    procedure :: bcast
end type

contains

subroutine mult_l(this,lat,res,work,alpha,beta)
    !this implementation without the transpose is only very slightly faster thatn the operator on this%H with transpose
    !so that in most cases it is not worth it to use this class if only mult_l is supposed to be optimized
    class(t_H_mkl_csr_mem),intent(in)   :: this
    type(lattice), intent(in)           :: lat
    real(8), intent(inout)              :: res(:)
    type(work_mode),intent(inout)       :: work
    real(8),intent(in),optional         :: alpha
    real(8),intent(in),optional         :: beta
    ! internal
    integer(C_int)                      :: stat
    real(8),pointer ,contiguous         :: modes(:)
    real(8)                             :: alp, bet
    integer                             :: work_size(N_work)

    if(present(alpha))then
        alp=alpha
    else
        alp=1.0d0
    endif
    if(present(beta))then
        bet=beta
    else
        bet=0.0d0
    endif
    Call this%mode_l%get_mode(lat,modes,work,work_size)
    stat=mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,alp,this%HT,this%descr,modes,bet,res)
#ifdef CPP_MKL
    if(stat/=SPARSE_STATUS_SUCCESS) STOP "failed MKL_SPBLAS routine in mult_l of m_H_sparse_mkl_mem"
#endif
    nullify(modes)
    work%offset=work%offset-work_size
end subroutine 



subroutine mult_l_disc(this,i_m,lat,N,ind_out,vec,ind_sum,ind_Mult,mat_mult,vec_mult)
    !Calculates the entries of the matrix * right vector product for the indices ind_out of the result vector
    USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
    ! input
    class(t_H_mkl_csr_mem), intent(in)  :: this
    type(lattice), intent(in)           :: lat
    integer, intent(in)                 :: i_m          !index of the comp's right mode in the inner dim_mode
    integer, intent(in)                 :: N            !number of indices to calculated
    integer, intent(in)                 :: ind_out(N)   !indices to be calculated
    ! output
    real(8),intent(out)                 :: vec(N) ! dim_modes_inner(this%mode_r%order(comp))
    !temporary data
    integer,intent(inout)               :: ind_sum (N+1)                !ind_mult index where the a vec-entry start and end
    integer,intent(inout)               :: ind_mult(N*this%col_max)     !indices of the left array which have non-vanishing contributions to get the ind entries of the vec/mat product
    real(8),intent(inout)               :: mat_mult(N*this%col_max)     !matrix entries corresponding to ind_mult
    real(8),intent(inout)               :: vec_mult(N*this%col_max)     !values of discontiguous left mode which has to be evaluated (indices of ind_mult)

    !some local indices/ loop variables
    integer ::  i,j, i_outer,ii

    !get matrix indices and values whose vec.mat product constitute the output vector
    ind_sum(1)=0
    ii=0
    do i=1,N
        i_outer=ind_out(i)
        do j=this%HT_outer(i_outer),this%HT_outer(i_outer+1)-1
            ii=ii+1
            mat_mult(ii)=this%HT_val(j)
            ind_mult(ii)=this%HT_inner(j)
        enddo
        ind_sum(i+1)=ii
    enddo

    !get the left vector values which are multiplied with the matrix
    Call this%mode_l%get_mode_disc(lat,ii,ind_mult(:ii),vec_mult(:ii))

    !multipy the matrix and left vector entries and sum together to the respective output entry
    vec_mult(:ii)=vec_mult(:ii)*mat_mult(:ii)
    do i=1,N
        vec(i)=sum(vec_mult(ind_sum(i)+1:ind_sum(i+1)))
    enddo
end subroutine 


subroutine optimize(this)
    class(t_H_mkl_csr_mem),intent(inout)   :: this

    Call this%t_H_mkl_csr%optimize()
    !THERE MIGHT BE SOMETHING TO DO HERE
    continue 
end subroutine

subroutine set_work_single(this,work,order)
    class(t_H_mkl_csr_mem),intent(inout)      :: this
    class(work_ham_single),intent(inout)    :: work 
    integer,intent(in)                      :: order
    integer     :: sizes(N_work)
    integer     :: dim_mode

    Call this%t_H_mkl_csr%set_work_single(work,order)
    if(this%col_max==0) ERROR STOP "cannot set work size of t_H_mkl_csr if col_max=0"
    Call this%get_work_size_single(sizes)
    Call work%set(sizes)
end subroutine

subroutine get_work_size_single(this,sizes)
    class(t_H_mkl_csr_mem),intent(in) :: this
    integer,intent(out)             :: sizes(N_work_single)
    integer                         :: sizes_n(N_work_single)

    Call this%t_H_mkl_csr%get_work_size_single(sizes_n)
    Call work_size_single(maxval(this%dim_r_single),this%col_max,sizes)
    sizes=max(sizes,sizes_n)
end subroutine

subroutine copy_child(this,Hout)
    class(t_H_mkl_csr_mem),intent(in)   :: this
    class(t_H_base),intent(inout)       :: Hout
    integer         ::  stat
    
    Call this%t_H_mkl_csr%copy_child(Hout)
    select type(Hout)
    class is(t_H_mkl_csr_mem)
        stat=mkl_sparse_copy(this%HT,this%descr,Hout%HT)
        if(stat/=SPARSE_STATUS_SUCCESS) STOP 'failed to copy Sparse_Matrix_T in m_H_sparse_mkl'
        Call set_auxiliaries(Hout)
    class default
        STOP "Cannot copy t_H_mkl_csr_mem type with Hamiltonian that is not a class of t_H_mkl_csr_mem"
    end select
end subroutine

subroutine add_child(this,H_in)
    class(t_H_mkl_csr_mem),intent(inout)  :: this
    class(t_H_base),intent(in)            :: H_in

    type(SPARSE_MATRIX_T)       :: tmp_H
    integer                     :: stat
    real(C_DOUBLE),parameter    :: alpha=1.0d0

    Call this%t_H_mkl_csr%add_child(H_in)
    select type(H_in)
    class is(t_H_mkl_csr_mem)
        tmp_H=this%HT
        stat=mkl_sparse_d_add(sparse_operation_non_transpose,tmp_H,alpha,h_in%HT,this%HT)
        if(stat/=SPARSE_STATUS_SUCCESS) STOP "add failed in mkl_inspector_csr"
        stat=mkl_sparse_destroy(tmp_H)
        if(stat/=SPARSE_STATUS_SUCCESS) STOP 'failed to destroy t_h_mkl_coo type in m_H_sparse_mkl'
        Call this%set_auxiliaries()
    class default
        STOP "Cannot add t_H_mkl_csr_mem type with Hamiltonian that is not a class of t_H_mkl_csr_mem"
    end select
end subroutine 

subroutine destroy_child(this)
    class(t_H_mkl_csr_mem),intent(inout)    :: this
    integer     ::  stat

    Call this%t_H_mkl_csr%destroy_child()
    if(this%is_set())then
        stat=mkl_sparse_destroy(this%HT)
        nullify(this%HT_val,this%HT_inner,this%HT_outer)
        this%col_max=0
        if(stat/=SPARSE_STATUS_SUCCESS) STOP 'failed to destroy t_h_mkl_csr type in m_H_sparse_mkl'
    endif
end subroutine

subroutine set_from_Hcoo(this,H_coo)
    class(t_H_mkl_csr_mem),intent(inout)  :: this
    type(t_H_coo),intent(inout)     :: H_coo

    type(t_H_coo)                   :: H_coo_copy

    type(SPARSE_MATRIX_T)   :: H
    integer(C_int)                  :: stat
    integer                         :: nnz
    real(C_DOUBLE),allocatable      :: val(:)
    integer(C_INT),allocatable      :: rowind(:),colind(:)


    !super wastefull to copy t_H_coo, but the fastest way to implement without changing the interface
    Call H_coo%copy(H_coo_copy)
    Call this%t_H_mkl_csr%set_from_Hcoo(H_coo_copy)

    !create mkl sparse matrix in coo format
    Call H_coo%pop_par(this%dimH,nnz,val,rowind,colind)
    stat=mkl_sparse_d_create_coo(H, SPARSE_INDEX_BASE_ONE , this%dimH(2) , this%dimH(1) , nnz , colind , rowind , val)
    if(stat/=SPARSE_STATUS_SUCCESS) STOP "failed to initialize MKL_SPBLAS matrix"

    !convert to mkl-sparse matrix in csr format
    stat = MKL_SPARSE_CONVERT_CSR(H,SPARSE_OPERATION_NON_TRANSPOSE,this%HT)
    if(stat /= 0) STOP "error setting H_ee sparse to CSR"

    !destroy coo-sparse matrix
    stat=MKL_SPARSE_DESTROY(H)
    if(stat /= 0) STOP "error destroying H"
   
    Call this%set_auxiliaries()
end subroutine 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!            MPI ROUTINES           !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine send(this,ithread,tag,com)
    use mpi_basic                
    class(t_H_mkl_csr_mem),intent(in) :: this
    integer,intent(in)              :: ithread
    integer,intent(in)              :: tag
    integer,intent(in)              :: com

#ifdef CPP_MPI
    Call this%t_H_mkl_csr%send(ithread,tag,com)
    STOP "IMPLEMENT ANALOGOUS TO t_H_mkl_csr"
#else
    continue
#endif
end subroutine

subroutine recv(this,ithread,tag,com)
    use mpi_basic                
    class(t_H_mkl_csr_mem),intent(inout)  :: this
    integer,intent(in)                  :: ithread
    integer,intent(in)                  :: tag
    integer,intent(in)                  :: com

#ifdef CPP_MPI
    Call this%t_H_mkl_csr%recv(ithread,tag,com)
    STOP "IMPLEMENT ANALOGOUS TO t_H_mkl_csr"
    Call this%set_auxiliaries()
#else
    continue
#endif
end subroutine

subroutine bcast(this,comm)
    use mpi_basic                
    class(t_H_mkl_csr_mem),intent(inout)  ::  this
    type(mpi_type),intent(in)           ::  comm
#ifdef CPP_MPI
    Call this%t_H_mkl_csr%bcast(comm)
    STOP "IMPLEMENT ANALOGOUS TO t_H_mkl_csr"
    if(.not.comm%ismas) Call this%set_auxiliaries()
#else
    continue
#endif
end subroutine 

subroutine distribute(this,comm)
    use mpi_basic                
    class(t_H_mkl_csr_mem),intent(inout)  ::  this
    type(mpi_type),intent(in)           ::  comm
#ifdef CPP_MPI
    Call this%t_H_mkl_csr%distribute(comm)
    STOP "IMPLEMENT ANALOGOUS TO t_H_mkl_csr"
    if(.not.comm%ismas) Call this%set_auxiliaries()
#else
    continue
#endif
end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           HELPER FUNCTIONS                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_auxiliaries(this)
    class(t_H_mkl_csr_mem),intent(inout)    :: this

    Call this%t_H_mkl_csr%set_auxiliaries()
    Call set_H_ptr(this)
    Call set_col_max(this)
end subroutine

subroutine set_H_ptr(this)
    class(t_H_mkl_csr_mem),intent(inout)    :: this

    type(C_PTR)     :: col,row,val
    integer(C_INT)  :: nnz

    nullify(this%HT_outer,this%HT_inner,this%HT_val)
    Call unpack_csr(this%HT,nnz,this%HT_outer,this%HT_inner,this%HT_val)
end subroutine

subroutine set_col_max(this)
    class(t_H_mkl_csr_mem),intent(inout)  :: this
    !variable to get multiplication for single evaluation (estimate sizes...)
    integer,allocatable               :: tmp(:)
    integer                           :: i

    allocate(tmp(size(this%HT_outer)-1))
    do i=1,size(tmp)
        tmp(i)=this%HT_outer(i+1)-this%HT_outer(i)
    enddo
    this%col_max=maxval(tmp)
end subroutine
#endif
end module m_H_mkl_csr_mem
