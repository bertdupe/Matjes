module m_H_sparse_mkl
#if defined(CPP_MKL_SPBLAS)
!Hamiltonian type specifications using MKL_SPARSE inspector mkl in csr 
!eval_single single energy evaluation is rather cumbersome...
use MKL_SPBLAS
use m_derived_types, only: lattice,number_different_order_parameters
use m_H_coo_based
use mkl_spblas_util, only: unpack_csr
USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_DOUBLE,C_INT


type,extends(t_H_coo_based) :: t_H_mkl_csr
    private
    !mkl parameters
    type(SPARSE_MATRIX_T)   :: H
    type(matrix_descr)      :: descr
    !pointers to Hamiltonian data handled by mkl
    real(C_DOUBLE),pointer  :: mat_val(:)
    integer(C_INT),pointer  :: mat_col(:),mat_row(:)
    integer                 :: nnz
    !helper variables
    integer                 :: col_max=0

contains
    !necessary t_H routines
    procedure :: eval_single

    procedure :: set_from_Hcoo

    procedure :: add_child 
    procedure :: destroy_child    
    procedure :: copy_child 

    procedure :: optimize
    procedure :: mult_r,mult_l
    procedure :: mult_l_cont,mult_r_cont
    procedure :: mult_l_disc,mult_r_disc

    !MPI
    procedure :: send
    procedure :: recv
    procedure :: bcast
    procedure :: distribute 
end type

interface t_H_mkl_csr
    procedure :: dummy_constructor
end interface 
 
private
public t_H,t_H_mkl_csr
contains 

type(t_H_mkl_csr) function dummy_constructor()
    !might want some initialization for H and descr, but should work without
    !continue 
end function 

subroutine mult_r(this,lat,res)
    !mult
    use m_derived_types, only: lattice
    class(t_H_mkl_csr),intent(in)   :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)   !result matrix-vector product
    ! internal
    integer(C_int)             :: stat
    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)
    real(C_DOUBLE),parameter   :: alpha=1.0d0,beta=0.0d0

    Call this%mode_r%get_mode(lat,modes,vec)
    if(size(res)/=this%dimH(1)) STOP "size of vec is wrong"
    res=0.0d0
    stat=mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,alpha,this%H,this%descr,modes,beta,res)
    if(stat/=SPARSE_STATUS_SUCCESS) STOP "failed MKL_SPBLAS routine in mult_r of m_H_sparse_mkl"
    if(allocated(vec)) deallocate(vec)
end subroutine 

subroutine mult_l(this,lat,res)
    use m_derived_types, only: lattice
    class(t_H_mkl_csr),intent(in)   :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)
    ! internal
    integer(C_int)             :: stat
    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)
    real(C_DOUBLE),parameter   :: alpha=1.0d0,beta=0.0d0

    Call this%mode_l%get_mode(lat,modes,vec)
    if(size(res)/=this%dimH(2)) STOP "size of vec is wrong"
    res=0.0d0
    stat=mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE,alpha,this%H,this%descr,modes,beta,res)
    if(stat/=SPARSE_STATUS_SUCCESS) STOP "failed MKL_SPBLAS routine in mult_l of m_H_sparse_mkl"
    if(allocated(vec)) deallocate(vec)
end subroutine 


subroutine mult_l_cont(this,bnd,vec,res)
    !multiply to the right with a continuous section of the right vector
    class(t_H_mkl_csr),intent(in)     :: this
    integer,intent(in)              :: bnd(2)
    real(8),intent(in)              :: vec(bnd(2)-bnd(1)+1)
    real(8),intent(inout)           :: res(:)   !result matrix-vector product

    ERROR STOP "IMPLEMENT"
end subroutine 

subroutine mult_r_cont(this,bnd,vec,res)
    !multiply to the right with a continuous section of the right vector
    class(t_H_mkl_csr),intent(in)     :: this
    integer,intent(in)              :: bnd(2)
    real(8),intent(in)              :: vec(bnd(2)-bnd(1)+1)
    real(8),intent(inout)           :: res(:)   !result matrix-vector product

    ERROR STOP "IMPLEMENT"
end subroutine 

subroutine mult_l_disc(this,N,ind,vec,res)
    !multiply to the right with a discontinuous section of the right vector
    class(t_H_mkl_csr),intent(in)     :: this
    integer,intent(in)              :: N
    integer,intent(in)              :: ind(N)
    real(8),intent(in)              :: vec(N)
    real(8),intent(inout)           :: res(:)   !result matrix-vector product

    ERROR STOP "IMPLEMENT"
end subroutine 

subroutine mult_r_disc(this,N,ind,vec,res)
    !multiply to the right with a discontinuous section of the right vector
    class(t_H_mkl_csr),intent(in)     :: this
    integer,intent(in)              :: N
    integer,intent(in)              :: ind(N)
    real(8),intent(in)              :: vec(N)
    real(8),intent(inout)           :: res(:)   !result matrix-vector product

    ERROR STOP "IMPLEMENT"
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
        Call this%copy_deriv(Hout)
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

        Call set_auxiliaries(this)
    class default
        STOP "Cannot add t_h_mkl_csr type with Hamiltonian that is not a class of t_h_mkl_csr"
    end select

end subroutine 

subroutine destroy_child(this)
    class(t_H_mkl_csr),intent(inout)    :: this
    integer     ::  stat

    if(this%is_set())then
        stat=mkl_sparse_destroy(this%H)
        nullify(this%mat_val,this%mat_col,this%mat_row)
        this%nnz=0
        this%col_max=0
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
   
    Call set_auxiliaries(this)
end subroutine 


subroutine eval_single(this,E,i_m,order,lat)
    use m_derived_types, only: lattice
    USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
    ! input
    class(t_H_mkl_csr), intent(in)  :: this
    type(lattice), intent(in)       :: lat
    integer, intent(in)             :: i_m
    integer, intent(in)             :: order
    ! output
    real(8), intent(out)            :: E

    integer,allocatable             :: ind(:)
    real(8),allocatable             :: vec(:)
    real(8),allocatable             :: vec_out(:)
    integer,allocatable             :: ind_out(:)
    real(8),allocatable             :: vec_l(:) !discontiguous array on left side
    integer                         :: N_out

    !some local indices/ loop variables
    integer ::  i,j, i_col
    integer :: ii

    Call this%mode_r%get_mode_single_disc(lat,1,i_m,ind,vec)
    N_out=size(ind)*this%col_max
    allocate(vec_out(N_out), ind_out(N_out))
    
    ii=0
    do i=1,size(ind)
        i_col=ind(i)
        do j=this%mat_row(i_col),this%mat_row(i_col+1)-1
            ii=ii+1
            vec_out(ii)=this%mat_val(j)*vec(i)
            ind_out(ii)=this%mat_col(j)
        enddo
    enddo

    Call this%mode_l%get_mode_disc(lat,ind_out(:ii),vec_l)
    E=DOT_PRODUCT(vec_l(:ii),vec_out(:ii))
end subroutine 

subroutine create_sparse_vec(i_m,modes,dim_mode,dim_H,vec)
    !returns a sparse matrix (vec) in csr-format which describes a sparse (dim_H,1) vector
    !with the ((i_m-1)*dim_mode+1:i_m*dim_mode)) values from modes
    type(SPARSE_MATRIX_T),intent(out) :: vec
    real(8),pointer,intent(in)        :: modes(:)
    integer,intent(in)                :: i_m
    integer,intent(in)                :: dim_mode,dim_H

    integer                 :: col(dim_mode)
    integer                 :: row(dim_mode)
    type(SPARSE_MATRIX_T)   :: vec_coo
    integer(C_int)          :: stat
    integer                 :: i

    col=1
    row=dim_mode*(i_m-1)
    do i=1,dim_mode
        row(i)=row(i)+i
    enddo

    stat=mkl_sparse_d_create_coo(vec_coo, SPARSE_INDEX_BASE_ONE , dim_H , 1 , dim_mode , row , col , modes((i_m-1)*dim_mode+1:i_m*dim_mode))
    if(stat/=SPARSE_STATUS_SUCCESS) ERROR STOP "mkl error"
    stat = MKL_SPARSE_CONVERT_CSR(vec_coo,SPARSE_OPERATION_NON_TRANSPOSE,vec)
    if(stat/=SPARSE_STATUS_SUCCESS) ERROR STOP "mkl error"
    stat=mkl_sparse_destroy(vec_coo)
    if(stat/=SPARSE_STATUS_SUCCESS) ERROR STOP "mkl error"
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

    Call this%set_deriv()

    Call set_auxiliaries(this)
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
        Call set_auxiliaries(this)
    endif
    nullify(acsr,ia,ja)
    if(.not.comm%ismas) Call this%set_deriv()
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
    Call set_auxiliaries(this)

    if(.not.comm%ismas) Call this%set_deriv()
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

    if(.not.comm%ismas) Call this%set_deriv()
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
    Call set_H_ptr(this)
    Call set_col_max(this)
end subroutine

subroutine set_H_ptr(this)
    class(t_H_mkl_csr),intent(inout)    :: this

    nullify(this%mat_row,this%mat_col,this%mat_val)
    Call unpack_csr(this%H,this%nnz,this%mat_row,this%mat_col,this%mat_val)
end subroutine

subroutine set_col_max(this)
    class(t_H_mkl_csr),intent(inout)    :: this
    !variable to get multiplication for single evaluation (estimate sizes...)
    integer,allocatable                 :: tmp(:)
    integer                             :: i

    allocate(tmp(size(this%mat_row)-1))
    do i=1,size(tmp)
        tmp(i)=this%mat_row(i+1)-this%mat_row(i)
    enddo
    this%col_max=maxval(tmp)
end subroutine


#endif
end module
