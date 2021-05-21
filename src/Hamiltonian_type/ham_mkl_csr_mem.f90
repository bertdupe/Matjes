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
    integer(C_INT),pointer  :: HT_col(:),HT_row(:)
    !helper variables
    integer                 :: col_max=0 !maximal number of entries per col

contains
    procedure :: set_from_Hcoo

    procedure :: add_child 
    procedure :: destroy_child    
    procedure :: copy_child 

    procedure :: optimize

    !might make sense to add full matrix multiplication to faster side with transpose

    procedure :: mult_l_single
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


subroutine mult_l_single(this,i_m,comp,lat,work,vec)
    !Calculates the entries of the left vector*matrix product which corresponds to the i_m's site of component comp of the right modes
    USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
    ! input
    class(t_H_mkl_csr_mem), intent(in)  :: this
    type(lattice), intent(in)           :: lat
    integer, intent(in)                 :: i_m           !index of the comp's right mode in the inner dim_mode
    integer, intent(in)                 :: comp          !component of right mode
    !temporary data
    type(work_ham_single),intent(inout) :: work          !data type containing the temporary data for this calculation to prevent constant allocations/deallocations
    ! output
    real(8),intent(inout)               :: vec(:)        !dim_modes_inner(this%mode_l%order(comp))

    integer :: i,j,i_row
    integer :: ii

#ifdef CPP_USE_WORK
    !temporary data slices
    integer,pointer,contiguous          :: ind_out(:)    !indices of all right mode entries which contain the order paramtere order of the site corresponding to i_m
    integer,pointer,contiguous          :: ind_sum(:)    !ind_mult index where the a vec-entry start and end
    integer,pointer,contiguous          :: ind_mult(:)   !indices of the left array which have non-vanishing contributions to get the ind entries of the vec/mat product
    real(8),pointer,contiguous          :: mat_mult(:)   !matrix entries corresponding to ind_mult
    real(8),pointer,contiguous          :: vec_l(:)      !values of discontiguous left mode which has to be evaluated (indices of ind_mult)


    !associate temporary arrays
    !!int vector slices
    ind_out (1:this%dim_r_single             )=>work%int_arr (1                               :this%dim_l_single                    )
    ind_sum (1:this%dim_r_single+1           )=>work%int_arr (1+this%dim_r_single             :this%dim_r_single* 2               +1)
    ind_mult(1:this%dim_r_single*this%col_max)=>work%int_arr (1+this%dim_r_single*2+1         :this%dim_r_single*(2+ this%col_max)+1)
    !!real vector slices
    vec_r   (1:this%dim_r_single*this%col_max)=>work%real_arr(1                               :this%dim_r_single*this%col_max                  )
    mat_mult(1:this%dim_r_single*this%col_max)=>work%real_arr(1+this%dim_r_single*this%col_max:this%dim_r_single*this%col_max*2                )
#else
    !temporary arrays
    integer                     :: ind_out (this%dim_l_single             )     !indices of all right mode entries which contain the order paramtere order of the site corresponding to i_m
    integer                     :: ind_sum (this%dim_l_single+1           )     !ind_mult index where the a vec-entry start and end
    integer                     :: ind_mult(this%dim_l_single*this%row_max)     !indices of the left array which have non-vanishing contributions to get the ind entries of the vec/mat product
    real(8)                     :: mat_mult(this%dim_l_single*this%row_max)     !matrix entries corresponding to ind_mult
    real(8)                     :: vec_l   (this%dim_l_single*this%row_max)     !values of discontiguous left mode which has to be evaluated (indices of ind_mult)
#endif
    !get indices of the output vector 
    Call this%mode_r%get_ind_site_expl(comp,i_m,this%dim_r_single,ind_out)

    !get matrix indices and values whose vec.mat product constitute the output vector
    ind_sum(1)=0
    ii=0
    do i=1,this%dim_r_single
        i_row=ind_out(i)
        do j=this%H_row(i_row),this%H_row(i_row+1)-1
            ii=ii+1
            mat_mult(ii)=this%H_val(j)
            ind_mult(ii)=this%H_col(j)
        enddo
        ind_sum(i+1)=ii
    enddo

    !get the right vector values which are multiplied with the matrix
    Call this%mode_l%get_mode_disc_expl(lat,ii,ind_mult(:ii),vec_l(:ii))

    !multipy the matrix and left vector entries and sum together to the respective output entry
    vec_l(:ii)=vec_l(:ii)*mat_mult(:ii)
    do i=1,this%dim_r_single
        vec(i)=sum(vec_l(ind_sum(i)+1:ind_sum(i+1)))
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
    integer     :: sizes(2)
    integer     :: dim_mode

    Call this%t_H_mkl_csr%set_work_single(work,order)
    if(this%col_max==0) ERROR STOP "cannot set work size of t_H_mkl_csr if col_max=0"
    if(.not.any(order==this%op_l))then
        if(any(order==this%op_r))then
            ERROR STOP "DECIDE HOW TO TREATE CASE WITHOUT ENTRY"
        else
            ERROR STOP "Hamiltonian has no component of considered single energy evaluation, take it out or consider it somehow else"
        endif
    endif
    Call this%get_work_size_single(sizes)
    Call work%set(sizes)
end subroutine

subroutine get_work_size_single(this,sizes)
    class(t_H_mkl_csr_mem),intent(in) :: this
    integer,intent(out)             :: sizes(N_work_single)
    integer                         :: sizes_n(N_work_single)

    Call this%t_H_mkl_csr%get_work_size_single(sizes_n)
    Call work_size_single(this%dim_l_single,this%row_max,sizes)
    Call work_size_single(this%dim_r_single,this%col_max,sizes)
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
        nullify(this%HT_val,this%HT_col,this%HT_row)
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

    nullify(this%HT_row,this%HT_col,this%HT_val)
    Call unpack_csr(this%HT,nnz,this%HT_row,this%HT_col,this%HT_val)
end subroutine

subroutine set_col_max(this)
    class(t_H_mkl_csr_mem),intent(inout)  :: this
    !variable to get multiplication for single evaluation (estimate sizes...)
    integer,allocatable               :: tmp(:)
    integer                           :: i

    allocate(tmp(size(this%HT_row)-1))
    do i=1,size(tmp)
        tmp(i)=this%HT_row(i+1)-this%HT_row(i)
    enddo
    this%col_max=maxval(tmp)
end subroutine
#endif
end module m_H_mkl_csr_mem
