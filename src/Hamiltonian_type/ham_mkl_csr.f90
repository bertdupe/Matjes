module m_H_mkl_csr
#if defined(CPP_MKL)
!Hamiltonian type specifications using MKL_SPARSE inspector mkl in csr 
use MKL_SPBLAS
use m_type_lattice, only: dim_modes_inner, lattice,number_different_order_parameters
use m_H_coo_based
use mkl_spblas_util, only: unpack_csr
use m_work_ham_single
USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_DOUBLE,C_INT
use,intrinsic :: ISO_FORTRAN_ENV, only: error_unit

private
public t_H,t_H_mkl_csr

type,extends(t_H_coo_based) :: t_H_mkl_csr
!    private
    !mkl parameters
    type(SPARSE_MATRIX_T)   :: H
    type(matrix_descr)      :: descr
    !pointers to Hamiltonian data handled by mkl (row major format)
    real(C_DOUBLE),pointer  :: H_val(:)
    integer(C_INT),pointer  :: H_inner(:),H_outer(:)
    integer                 :: nnz
contains
    !necessary t_H routines

    procedure :: set_from_Hcoo

    procedure :: add_child 
    procedure :: destroy_child    
    procedure :: copy_child 

    procedure :: optimize
    procedure :: mult_r,mult_l

    procedure :: mult_r_disc, mult_l_disc

    procedure :: set_work_single
    procedure :: get_work_size_single

    !utility
    procedure :: set_auxiliaries

    !MPI
    procedure :: send
    procedure :: recv
    procedure :: bcast
    procedure :: distribute 
end type

interface t_H_mkl_csr
    procedure :: dummy_constructor
end interface 
 
contains 

subroutine set_work_single(this,work,order)
    class(t_H_mkl_csr),intent(inout)        :: this
    class(work_ham_single),intent(inout)    :: work 
    integer,intent(in)                      :: order
    integer     :: sizes(2)
    integer     :: dim_mode

    if(.not.this%is_set()) ERROR STOP "cannot set work size of hamiltonian if it is not set"
    if(this%row_max==0) ERROR STOP "cannot set work size of t_H_mkl_csr if row_max==0"
    Call this%get_work_size_single(sizes)
    Call work%set(sizes)
end subroutine

subroutine get_work_size_single(this,sizes)
    class(t_H_mkl_csr),intent(in)   :: this
    integer,intent(out)             :: sizes(2)

    Call work_size_single(maxval(this%dim_l_single),this%row_max,sizes)
end subroutine

type(t_H_mkl_csr) function dummy_constructor()
    !might want some initialization for H and descr, but should work without
    !continue 
end function 

subroutine mult_r(this,lat,res,work,alpha,beta)
    class(t_H_mkl_csr),intent(in)   :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)   !result matrix-vector product
    type(work_mode),intent(inout)   :: work
    real(8),intent(in),optional     :: alpha
    real(8),intent(in),optional     :: beta
    ! internal
    integer(C_int)                  :: stat
    real(8),pointer ,contiguous     :: modes(:)
    real(8)                         :: alp, bet
    integer                         :: work_size(N_work)

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

    Call this%mode_r%get_mode(lat,modes,work,work_size)
    stat=mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,alp,this%H,this%descr,modes,bet,res)
#ifdef CPP_DEBUG
    if(stat/=SPARSE_STATUS_SUCCESS) STOP "failed MKL_SPBLAS routine in mult_r of m_H_sparse_mkl"
#endif
    nullify(modes)
    work%offset=work%offset-work_size
end subroutine 

subroutine mult_l(this,lat,res,work,alpha,beta)
    class(t_H_mkl_csr),intent(in)   :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)
    type(work_mode),intent(inout)   :: work
    real(8),intent(in),optional     :: alpha
    real(8),intent(in),optional     :: beta
    ! internal
    integer(C_int)                  :: stat
    real(8),pointer ,contiguous     :: modes(:)
    real(8)                         :: alp, bet
    integer                         :: work_size(N_work)

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
    stat=mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE,alp,this%H,this%descr,modes,bet,res)
#ifdef CPP_DEBUG
    if(stat/=SPARSE_STATUS_SUCCESS) STOP "failed MKL_SPBLAS routine in mult_l of m_H_sparse_mkl"
#endif
    nullify(modes)
    work%offset=work%offset-work_size
end subroutine 

subroutine optimize(this)
    class(t_H_mkl_csr),intent(inout)   :: this
    integer     :: stat

    stat=mkl_sparse_order(this%H)
    stat=mkl_sparse_set_dotmv_hint ( this%H , SPARSE_OPERATION_NON_TRANSPOSE , this%descr , 100000 )
    stat=mkl_sparse_optimize(this%H)
end subroutine

subroutine copy_child(this,Hout)
    class(t_H_mkl_csr),intent(in)   :: this
    class(t_H_base),intent(inout)   :: Hout
    integer         ::  stat
    
    select type(Hout)
    class is(t_H_mkl_csr)
        stat=mkl_sparse_copy(this%H,this%descr,Hout%H )
        Hout%descr=this%descr
        if(stat/=SPARSE_STATUS_SUCCESS) STOP 'failed to copy Sparse_Matrix_T in m_H_sparse_mkl'
        Call set_auxiliaries(Hout)
    class default
        STOP "Cannot copy t_h_mkl_csr type with Hamiltonian that is not a class of t_h_mkl_csr"
    end select
end subroutine

subroutine add_child(this,H_in)
    class(t_H_mkl_csr),intent(inout)    :: this
    class(t_H_base),intent(in)          :: H_in
    
    type(SPARSE_MATRIX_T)       :: tmp_H
    integer                     :: stat
    real(C_DOUBLE),parameter    :: alpha=1.0d0


    select type(H_in)
    class is(t_H_mkl_csr)
        tmp_H=this%H
        stat=mkl_sparse_d_add(sparse_operation_non_transpose,tmp_H,alpha,h_in%h,this%H)
        if(stat/=SPARSE_STATUS_SUCCESS) STOP "add failed in mkl_inspector_csr"
        stat=mkl_sparse_destroy(tmp_H)
        if(stat/=SPARSE_STATUS_SUCCESS) STOP 'failed to destroy t_h_mkl_coo type in m_H_sparse_mkl'

        Call this%set_auxiliaries()
    class default
        STOP "Cannot add t_h_mkl_csr type with Hamiltonian that is not a class of t_h_mkl_csr"
    end select

end subroutine 

subroutine destroy_child(this)
    class(t_H_mkl_csr),intent(inout)    :: this
    integer     ::  stat

    if(this%is_set())then
        stat=mkl_sparse_destroy(this%H)
        nullify(this%H_val,this%H_inner,this%H_outer)
        this%nnz=0
        this%row_max=0
        if(stat/=SPARSE_STATUS_SUCCESS) STOP 'failed to destroy t_h_mkl_csr type in m_H_sparse_mkl'
    endif
end subroutine

subroutine set_from_Hcoo(this,H_coo)
    class(t_H_mkl_csr),intent(inout)    :: this
    type(t_H_coo),intent(inout)         :: H_coo

    !local
    integer                 :: nnz
    integer(C_int)          :: stat
    type(SPARSE_MATRIX_T)   :: H
    real(C_DOUBLE),allocatable     :: val(:)
    integer,allocatable     :: rowind(:),colind(:)

    !set descriptions
    this%descr%type=SPARSE_MATRIX_TYPE_GENERAL 
    this%descr%diag=SPARSE_DIAG_NON_UNIT
    this%descr%mode=SPARSE_FILL_MODE_LOWER

    !create mkl sparse matrix in coo format
    Call H_coo%pop_par(this%dimH,nnz,val,rowind,colind)
    stat=mkl_sparse_d_create_coo(H, SPARSE_INDEX_BASE_ONE , this%dimH(1) , this%dimH(2) , nnz , rowind , colind , val)
    if(stat/=SPARSE_STATUS_SUCCESS) STOP "failed to initialize MKL_SPBLAS matrix"

    !convert to mkl-sparse matrix in csr format
    stat = MKL_SPARSE_CONVERT_CSR(H,SPARSE_OPERATION_NON_TRANSPOSE,this%H)
    if(stat /= 0) STOP "error setting H_ee sparse to CSR"

    !destroy coo-sparse matrix
    stat=MKL_SPARSE_DESTROY(H)
    if(stat /= 0) STOP "error destroying H"
   
    Call this%set_auxiliaries()
end subroutine 

!subroutine eval_single(this,E,i_m,order,lat,work)
!    USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
!    ! input
!    class(t_H_mkl_csr), intent(in)      :: this
!    type(lattice), intent(in)           :: lat
!    integer, intent(in)                 :: i_m
!    integer, intent(in)                 :: order
!    ! output
!    real(8), intent(out)                :: E
!    !temporary data
!    type(work_ham_single),intent(inout) ::  work    !data type containing the temporary data for this calculation to prevent constant allocations/deallocations
!    !temporary data slices
!    integer,pointer,contiguous          :: ind(:)       !indices of all left mode entries which contain the order paramtere order of the site corresponding to i_m
!    real(8),pointer,contiguous          :: vec(:)       !values corresponding to ind
!    integer,pointer,contiguous          :: ind_out(:)   !indices of the result array multipling the vector (ind/vec) to the matrix
!    real(8),pointer,contiguous          :: vec_out(:)   !values corresponding to ind_out
!    real(8),pointer,contiguous          :: vec_mult(:)     !values of discontiguous mode array on right side (indices of ind_out)
!
!    !some local indices/ loop variables
!    integer ::  i,j, i_row
!    integer :: ii
!#define _dim_ this%dim_l_single(order)
!
!    !associate temporary arrays
!    ind     (1:_dim_             )=>work%int_arr (1                       :_dim_                   )
!    ind_out (1:_dim_*this%row_max)=>work%int_arr (1+_dim_                 :_dim_*(1+  this%row_max))
!    vec     (1:_dim_             )=>work%real_arr(1                       :_dim_                   )
!    vec_out (1:_dim_*this%row_max)=>work%real_arr(1+_dim_                 :_dim_*(1+  this%row_max))
!    vec_mult(1:_dim_*this%row_max)=>work%real_arr(1+_dim_*(1+this%row_max):_dim_*(1+2*this%row_max))
!
!    !get left mode corresponding to site i_m of order order
!    Call this%mode_l%get_mode_single(lat,1,i_m,_dim_,ind,vec)    !get this to work with different orders (1 is not order here but component of left mode)
!
!    !Calculate left mode,matrix product only for the necessary discontiguous mode-indices
!    ii=0
!    do i=1,_dim_
!        i_row=ind(i)
!        do j=this%H_outer(i_row),this%H_outer(i_row+1)-1
!            ii=ii+1
!            vec_out(ii)=this%H_val(j)*vec(i)
!            ind_out(ii)=this%H_inner(j)
!        enddo
!    enddo
!    
!    !get right mode for indices of the vec/mat product 
!    Call this%mode_r%get_mode_disc(lat,ii,ind_out(:ii),vec_mult(:ii))
!
!    !Get the energy
!    E=DOT_PRODUCT(vec_out(:ii),vec_mult(:ii))
!    nullify(ind,ind_out,vec,vec_out,vec_mult)
!#undef _dim_
!end subroutine 

subroutine mult_r_disc(this,i_m,lat,N,ind_out,vec,ind_sum,ind_Mult,mat_mult,vec_mult)
    !Calculates the entries of the left vector * matrix product for the indices ind_out of the result vector
    ! input
    class(t_H_mkl_csr), intent(in)      :: this
    type(lattice), intent(in)           :: lat
    integer, intent(in)                 :: i_m          !index of the comp's right mode in the inner dim_mode
    integer, intent(in)                 :: N            !number of indices to calculated
    integer, intent(in)                 :: ind_out(N)   !indices to be calculated
    ! output
    real(8),intent(out)                 :: vec(N) ! dim_modes_inner(this%mode_r%order(comp))
    !temporary data
    integer,intent(inout)               :: ind_sum (N+1)                !ind_mult index where the a vec-entry start and end
    integer,intent(inout)               :: ind_mult(N*this%row_max)     !indices of the left array which have non-vanishing contributions to get the ind entries of the vec/mat product
    real(8),intent(inout)               :: mat_mult(N*this%row_max)     !matrix entries corresponding to ind_mult
    real(8),intent(inout)               :: vec_mult(N*this%row_max)     !values of discontiguous right mode which has to be evaluated (indices of ind_mult)

    !some local indices/ loop variables
    integer ::  i,j, i_outer,ii

    !get matrix indices and values whose vec.mat product constitute the output vector
    ind_sum(1)=0
    ii=0
    do i=1,N
        i_outer=ind_out(i)
        do j=this%H_outer(i_outer),this%H_outer(i_outer+1)-1
            ii=ii+1
            mat_mult(ii)=this%H_val(j)
            ind_mult(ii)=this%H_inner(j)
        enddo
        ind_sum(i+1)=ii
    enddo

    !get the right vector values which are multiplied with the matrix
    Call this%mode_r%get_mode_disc(lat,ii,ind_mult(:ii),vec_mult(:ii))

    !multipy the matrix and right vector entries and sum together to the respective output entry
    vec_mult(:ii)=vec_mult(:ii)*mat_mult(:ii)
    do i=1,N
        vec(i)=sum(vec_mult(ind_sum(i)+1:ind_sum(i+1)))
    enddo
end subroutine 

subroutine mult_l_disc(this,i_m,lat,N,ind_out,vec,ind_sum,ind_Mult,mat_mult,vec_mult)
    !Calculates the entries of the matrix * right vector product for the indices ind_out of the result vector
    ! input
    class(t_H_mkl_csr), intent(in)      :: this
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
    real(8),intent(inout)               :: vec_mult(N*this%col_max)     !values of discontiguous right mode which has to be evaluated (indices of ind_mult)

    write(error_unit,'(//A)') "Trying to call mult_l_dist from type t_H_mkl_csr."
    write(error_unit,'(A)')   "This can not be done efficiently without the transpose as implemented in t_H_mkl_csr_mem."
    write(error_unit,'(A)')   "Please choose the Hamiltonian implementation including the transpose (Hamiltonian_mode)."
    ERROR STOP
end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!            MPI ROUTINES           !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine send(this,ithread,tag,com)
    use mpi_basic                
    class(t_H_mkl_csr),intent(in)   :: this
    integer,intent(in)              :: ithread
    integer,intent(in)              :: tag
    integer,intent(in)              :: com
#ifdef CPP_MPI
    integer(C_int),pointer          :: ia(:),ja(:)
    real(C_DOUBLE),pointer          :: val(:)
    
    integer     :: nnz
    integer     :: ierr

    Call this%send_base(ithread,tag,com)

    Call unpack_csr(this%H,nnz,ia,ja,val)
    Call MPI_Send(nnz, 1, MPI_INT, ithread, tag,  com,  ierr)
    Call MPI_Send(ia , this%dimH(1)+1, MPI_INT,              ithread, tag,  com, ierr)
    Call MPI_Send(ja , nnz,            MPI_INT,              ithread, tag,  com, ierr)
    Call MPI_Send(val, nnz,            MPI_DOUBLE_PRECISION, ithread, tag,  com, ierr)

    nullify(ia,ja,val)
#else
    continue
#endif
end subroutine

subroutine recv(this,ithread,tag,com)
    use mpi_basic                
    class(t_H_mkl_csr),intent(inout)   :: this
    integer,intent(in)              :: ithread
    integer,intent(in)              :: tag
    integer,intent(in)              :: com

#ifdef CPP_MPI

    integer(C_int),allocatable     :: ia(:),ja(:)
    real(C_DOUBLE),allocatable     :: val(:)
    
    integer     :: nnz
    integer     :: i,ierr
    type(SPARSE_MATRIX_T) :: H_local
    type(matrix_descr)    :: descr
    integer     :: stat(MPI_STATUS_SIZE)

    Call this%recv_base(ithread,tag,com)

    Call MPI_Recv(nnz, 1, MPI_INT, ithread, tag,  com, stat, ierr)
    allocate(ia(this%dimH(1)+1),ja(nnz),val(nnz))
    Call MPI_Recv(ia , this%dimH(1)+1, MPI_INT,              ithread, tag,  com, stat, ierr)
    Call MPI_Recv(ja , nnz,            MPI_INT,              ithread, tag,  com, stat, ierr)
    Call MPI_Recv(val, nnz,            MPI_DOUBLE_PRECISION, ithread, tag,  com, stat, ierr)
    ierr=mkl_sparse_d_create_csr(H_local, SPARSE_INDEX_BASE_ONE , this%dimH(1) , this%dimH(2), ia(1:size(ia)-1), ia(2:size(ia)), ja, val)
    if(ierr/=SPARSE_STATUS_SUCCESS) ERROR STOP 'failed to create local mkl sparse matrix'
    descr%type=SPARSE_MATRIX_TYPE_GENERAL 
    descr%diag=SPARSE_DIAG_NON_UNIT
    descr%mode=SPARSE_FILL_MODE_LOWER

    ierr= mkl_sparse_copy ( H_local, descr , this%H)
    if(ierr/=SPARSE_STATUS_SUCCESS) ERROR STOP 'failed to copy mkl sparse Hamiltonian'
    this%descr%type=SPARSE_MATRIX_TYPE_GENERAL 
    this%descr%diag=SPARSE_DIAG_NON_UNIT
    this%descr%mode=SPARSE_FILL_MODE_LOWER
    Call this%optimize()

    Call this%set_auxiliaries()
#else
    continue
#endif
end subroutine

subroutine bcast(this,comm)
    use mpi_basic
    use mpi_util,only: bcast_util => bcast
    class(t_H_mkl_csr),intent(inout)    ::  this
    type(mpi_type),intent(in)           ::  comm
#ifdef CPP_MPI
    real(C_DOUBLE),pointer      :: acsr(:)
    integer(C_INT),pointer      :: ia(:),ja(:)
    integer                     :: nnz
    integer                     :: ierr
    type(SPARSE_MATRIX_T)       :: H_tmp


    Call this%bcast_base(comm)
    nullify(acsr,ia,ja)
    if(comm%ismas)then
        Call unpack_csr(this%H,nnz,ia,ja,acsr) 
    endif
    Call MPI_Bcast(nnz, 1, MPI_INTEGER, comm%mas, comm%com,ierr)
    if(.not.comm%ismas)then
        allocate(acsr(nnz),ja(nnz),ia(this%dimH(1)+1)) 
    endif
    Call bcast_util(ia,comm)
    Call bcast_util(ja,comm)
    Call bcast_util(acsr,comm)
    Call bcast_util(this%descr%type,comm)
    Call bcast_util(this%descr%mode,comm)
    Call bcast_util(this%descr%diag,comm)
    if(.not.comm%ismas)then
        ierr = mkl_sparse_d_create_csr( H_tmp, SPARSE_INDEX_BASE_ONE,this%dimH(1), this%dimH(2), ia(1:size(ia)-1), ia(2:size(ia)),ja ,acsr)
        if(ierr /= 0) ERROR STOP "FAILED TO CREATE CHILD MKL SPARSE HAMILTONIAN"
        !copy to savely deallocate data arrays (one could just leave it allocated since the memory isn't really lost, but it feels weird)
        ierr = mkl_sparse_copy(H_tmp, this%descr, this%H) 
        Call this%optimize()
        deallocate(acsr,ja,ia)
    endif
    nullify(acsr,ia,ja)
    if(.not.comm%ismas) Call this%set_auxiliaries()
#else
    continue
#endif
end subroutine 

#if 1
subroutine distribute(this,comm)
    use mpi_basic                
    use mpi_util!,only: bcast_util => bcast
    class(t_H_mkl_csr),intent(inout)        ::  this
    type(mpi_type),intent(in)       ::  comm
    real(C_DOUBLE),pointer          :: val_base(:)
    integer(C_INT),pointer          :: ia_base(:),ja_base(:)
#ifdef CPP_MPI
    integer                         :: nnz_base
    integer     ::  cnt(comm%Np),displ(comm%Np)

    integer(C_int),allocatable     :: ia(:),ja(:)
    real(C_DOUBLE),allocatable     :: val(:)

    integer(C_INT),target  :: tmpi(1)
    real(C_DOUBLE),target  :: tmpr(1)

    integer     ::  i,ierr
    type(SPARSE_MATRIX_T) :: H_local
    type(matrix_descr)    :: descr

    Call this%bcast_base(comm)
    if(comm%ismas)then
        Call unpack_csr(this%H,nnz_base,ia_base,ja_base,val_base)
    else
        val_base=> tmpr; ja_base=> tmpi
    endif
    Call bcast(nnz_base,comm)
    cnt=nnz_base/comm%Np
    forall(i=1:modulo(nnz_base,comm%Np)) cnt(i)=cnt(i)+1
    displ=[(sum(cnt(:i-1)),i=1,comm%np)]

    allocate(ia(this%dimH(1)+1),ja(cnt(comm%id+1)),val(cnt(comm%id+1)))
    if(comm%ismas) ia=ia_base
    Call bcast(ia,comm)
    ia=ia-displ(comm%id+1)
    ia=max(ia,1)
    ia=min(ia,cnt(comm%id+1)+1)

    if(comm%ismas)then
        Call MPI_Scatterv(ja_base,  cnt, displ, MPI_INT,              ja,  cnt(comm%id+1), MPI_INT,    comm%mas, comm%com, ierr)
        Call MPI_Scatterv(val_base, cnt, displ, MPI_DOUBLE_PRECISION, val, cnt(comm%id+1), MPI_DOUBLE, comm%mas, comm%com, ierr)
    else
        Call MPI_Scatterv(ja_base,  cnt, displ, MPI_INT,              ja,  cnt(comm%id+1), MPI_INT,    comm%mas, comm%com, ierr)
        Call MPI_Scatterv(val_base, cnt, displ, MPI_DOUBLE_PRECISION, val, cnt(comm%id+1), MPI_DOUBLE, comm%mas, comm%com, ierr)
    endif
    if(comm%ismas)then
        ierr=mkl_sparse_destroy(this%H)
        if(ierr/=SPARSE_STATUS_SUCCESS) ERROR STOP 'failed to destroy t_h_mkl_csr type in m_H_sparse_mkl'
    endif
    ierr=mkl_sparse_d_create_csr(H_local, SPARSE_INDEX_BASE_ONE , this%dimH(1) , this%dimH(2), ia(1:size(ia)-1), ia(2:size(ia)), ja, val)
    if(ierr/=SPARSE_STATUS_SUCCESS) ERROR STOP 'failed to create local mkl sparse matrix'
    descr%type=SPARSE_MATRIX_TYPE_GENERAL 
    descr%diag=SPARSE_DIAG_NON_UNIT
    descr%mode=SPARSE_FILL_MODE_LOWER

    ierr= mkl_sparse_copy ( H_local, descr , this%H)
    if(ierr/=SPARSE_STATUS_SUCCESS) ERROR STOP 'failed to copy mkl sparse Hamiltonian'
    this%descr%type=SPARSE_MATRIX_TYPE_GENERAL 
    this%descr%diag=SPARSE_DIAG_NON_UNIT
    this%descr%mode=SPARSE_FILL_MODE_LOWER
    Call this%optimize()

    if(.not.comm%ismas) Call this%set_auxiliaries()
#else
    continue
#endif
end subroutine 
#else
subroutine distribute(this,comm)
    use mpi_basic                
    use mpi_util!,only: bcast_util => bcast
    class(t_H_mkl_csr),intent(inout)        ::  this
    type(mpi_type),intent(in)       ::  comm
    real(C_DOUBLE),pointer          :: val_base(:)
    integer(C_INT),pointer          :: ia_base(:),ja_base(:)
#ifdef CPP_MPI
    integer                         :: nnz_base
    integer     ::  cnt(comm%Np),displ(comm%Np)

    integer(C_int),allocatable     :: ia(:),ja(:)
    real(C_DOUBLE),allocatable     :: val(:)

    integer(C_INT),target  :: tmpi(1)
    real(C_DOUBLE),target  :: tmpr(1)

    integer     ::  i,ierr
    type(SPARSE_MATRIX_T) :: H_local
    type(matrix_descr)    :: descr

    Call this%bcast_base(comm)
    if(comm%ismas)then
        Call unpack_csr(this%H,nnz_base,ia_base,ja_base,val_base)
    else
        val_base=> tmpr; ja_base=> tmpi
    endif
    Call bcast(nnz_base,comm)
    cnt=nnz_base/comm%Np
    forall(i=1:modulo(nnz_base,comm%Np)) cnt(i)=cnt(i)+1
    displ=[(sum(cnt(:i-1)),i=1,comm%np)]

    allocate(ia(this%dimH(1)+1),ja(cnt(comm%id+1)),val(cnt(comm%id+1)))
    if(comm%ismas) ia=ia_base
    Call bcast(ia,comm)
    ia=ia-displ(comm%id+1)
    ia=max(ia,1)
    ia=min(ia,cnt(comm%id+1)+1)

    if(comm%ismas)then
        Call MPI_Scatterv(ja_base,  cnt, displ, MPI_INT,              ja,  cnt(comm%id+1), MPI_INT,    comm%mas, comm%com, ierr)
        Call MPI_Scatterv(val_base, cnt, displ, MPI_DOUBLE_PRECISION, val, cnt(comm%id+1), MPI_DOUBLE, comm%mas, comm%com, ierr)
    else
        Call MPI_Scatterv(ja_base,  cnt, displ, MPI_INT,              ja,  cnt(comm%id+1), MPI_INT,    comm%mas, comm%com, ierr)
        Call MPI_Scatterv(val_base, cnt, displ, MPI_DOUBLE_PRECISION, val, cnt(comm%id+1), MPI_DOUBLE, comm%mas, comm%com, ierr)
    endif
    if(comm%ismas)then
        ierr=mkl_sparse_destroy(this%H)
        if(ierr/=SPARSE_STATUS_SUCCESS) ERROR STOP 'failed to destroy t_h_mkl_csr type in m_H_sparse_mkl'
    endif
    ierr=mkl_sparse_d_create_csr(H_local, SPARSE_INDEX_BASE_ONE , this%dimH(1) , this%dimH(2), ia(1:size(ia)-1), ia(2:size(ia)), ja, val)
    if(ierr/=SPARSE_STATUS_SUCCESS) ERROR STOP 'failed to create local mkl sparse matrix'
    descr%type=SPARSE_MATRIX_TYPE_GENERAL 
    descr%diag=SPARSE_DIAG_NON_UNIT
    descr%mode=SPARSE_FILL_MODE_LOWER

    ierr= mkl_sparse_copy ( H_local, descr , this%H)
    if(ierr/=SPARSE_STATUS_SUCCESS) ERROR STOP 'failed to copy mkl sparse Hamiltonian'
    this%descr%type=SPARSE_MATRIX_TYPE_GENERAL 
    this%descr%diag=SPARSE_DIAG_NON_UNIT
    this%descr%mode=SPARSE_FILL_MODE_LOWER
    Call this%optimize()

    if(.not.comm%ismas) call this%set_auxiliaries()
#else
    continue
#endif
end subroutine 
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           HELPER FUNCTIONS                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_auxiliaries(this)
    class(t_H_mkl_csr),intent(inout)    :: this

!    Call this%set_deriv()  !probably not necessary here
    Call set_H_ptr(this)
    Call set_row_max(this)
end subroutine

subroutine set_H_ptr(this)
    class(t_H_mkl_csr),intent(inout)    :: this

    nullify(this%H_outer,this%H_inner,this%H_val)
    Call unpack_csr(this%H,this%nnz,this%H_outer,this%H_inner,this%H_val)
end subroutine

subroutine set_row_max(this)
    class(t_H_mkl_csr),intent(inout)    :: this
    !variable to get multiplication for single evaluation (estimate sizes...)
    integer,allocatable                 :: tmp(:)
    integer                             :: i

    allocate(tmp(size(this%H_outer)-1))
    do i=1,size(tmp)
        tmp(i)=this%H_outer(i+1)-this%H_outer(i)
    enddo
    this%row_max=maxval(tmp)
end subroutine

#endif
end module
