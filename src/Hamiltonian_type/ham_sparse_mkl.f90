
module m_H_sparse_mkl
#if defined(CPP_MKL_SPBLAS)
!Hamiltonian type specifications using MKL_SPARSE inspector mkl in csr 
!eval_single single energy evaluation is rather cumbersome...
use m_H_type
use MKL_SPBLAS
use m_H_coo
use m_derived_types, only: lattice
USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_DOUBLE,C_INT


type,extends(t_H) :: t_H_mkl_csr
    private
    type(SPARSE_MATRIX_T) :: H
    type(matrix_descr)    :: descr
contains
    !necessary t_H routines
    procedure :: eval_single
    procedure :: init      
    procedure :: init_1    
    procedure :: init_mult_2   

    procedure :: add_child 
    procedure :: destroy_child    
    procedure :: copy_child 

    procedure :: optimize
    procedure :: mult_r,mult_l
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

    Call lat%point_order(this%op_r,this%dimH(2),modes,vec)

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

    Call lat%point_order(this%op_l,this%dimH(1),modes,vec)

    if(size(res)/=this%dimH(2)) STOP "size of vec is wrong"
    res=0.0d0
    stat=mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE,alpha,this%H,this%descr,modes,beta,res)
    if(stat/=SPARSE_STATUS_SUCCESS) STOP "failed MKL_SPBLAS routine in mult_l of m_H_sparse_mkl"
    if(allocated(vec)) deallocate(vec)
    
end subroutine 


subroutine optimize(this)
    class(t_H_mkl_csr),intent(inout)   :: this

    integer     :: stat

    stat=mkl_sparse_set_dotmv_hint ( this%H , SPARSE_OPERATION_NON_TRANSPOSE , this%descr , 100000 )
    stat=mkl_sparse_optimize(this%H)
    
end subroutine

subroutine init_1(this,line,Hval,Hval_ind,order,lat)
    use m_derived_types, only: lattice
    class(t_H_mkl_csr),intent(inout)    :: this

    type(lattice),intent(in)        :: lat
    integer,intent(in)              :: order(2)
    real(8),intent(in)              :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
    integer,intent(in)              :: Hval_ind(:,:)
    integer,intent(in)              :: line(:,:)

    !local
    type(t_H_coo)           :: H_coo

    if(this%is_set()) STOP "cannot set hamiltonian as it is already set"
    Call H_coo%init_1(line,Hval,Hval_ind,order,lat)
    Call set_from_Hcoo(this,H_coo,lat)
end subroutine 


subroutine init_mult_2(this,connect,Hval,Hval_ind,op_l,op_r,lat)
    use m_derived_types, only: lattice
    class(t_H_mkl_csr),intent(inout)    :: this

    type(lattice),intent(in)        :: lat
    integer,intent(in)              :: op_l(:),op_r(:)
    real(8),intent(in)              :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
    integer,intent(in)              :: Hval_ind(:,:)
    integer,intent(in)              :: connect(:,:)

    !local
    type(t_H_coo)           :: H_coo

    if(this%is_set()) STOP "cannot set hamiltonian as it is already set"
    Call H_coo%init_mult_2(connect,Hval,Hval_ind,op_l,op_r,lat)
    Call set_from_Hcoo(this,H_coo,lat)
end subroutine 


subroutine copy_child(this,Hout)
    class(t_H_mkl_csr),intent(in)   :: this
    class(t_H),intent(inout)        :: Hout
    integer         ::  stat
    
    select type(Hout)
    class is(t_H_mkl_csr)
        stat=mkl_sparse_copy(this%H,this%descr,Hout%H )
        Hout%descr=this%descr
        if(stat/=SPARSE_STATUS_SUCCESS) STOP 'failed to copy Sparse_Matrix_T in m_H_sparse_mkl'
    class default
        STOP "Cannot copy t_h_mkl_csr type with Hamiltonian that is not a class of t_h_mkl_csr"
    end select
end subroutine

subroutine add_child(this,H_in)
    class(t_H_mkl_csr),intent(inout)    :: this
    class(t_H),intent(in)               :: H_in
    
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
    class default
        STOP "Cannot add t_h_mkl_csr type with Hamiltonian that is not a class of t_h_mkl_csr"
    end select

end subroutine 

subroutine destroy_child(this)
    class(t_H_mkl_csr),intent(inout)    :: this
    integer     ::  stat

    if(this%is_set())then
        stat=mkl_sparse_destroy(this%H)
        if(stat/=SPARSE_STATUS_SUCCESS) STOP 'failed to destroy t_h_mkl_csr type in m_H_sparse_mkl'
    endif
end subroutine

subroutine init(this,energy_in,lat)
    use m_derived_types, only: operator_real_order_N,lattice
    class(t_H_mkl_csr),intent(inout)    :: this
    type(operator_real_order_N)         :: energy_in
    type(lattice),intent(in)            :: lat

    !local
    type(t_H_coo)           :: H_coo

    if(this%is_set()) STOP "cannot set hamiltonian as it is already set"
    Call H_coo%init(energy_in,lat)
    Call set_from_Hcoo(this,H_coo,lat)

end subroutine 

subroutine set_from_Hcoo(this,H_coo,lat)
    type(t_H_coo),intent(inout)         :: H_coo
    type(t_H_mkl_csr),intent(inout)     :: this
    type(lattice),intent(in)            :: lat

    !local
    integer                 :: nnz
    integer(C_int)          :: stat
    type(SPARSE_MATRIX_T)   :: H
    real(C_DOUBLE),allocatable     :: val(:)
    integer,allocatable     :: rowind(:),colind(:)

    Call H_coo%pop_par(this%dimH,nnz,val,rowind,colind)
    stat=mkl_sparse_d_create_coo(H, SPARSE_INDEX_BASE_ONE , this%dimH(1) , this%dimH(2) , nnz , rowind , colind , val)

    if(stat/=SPARSE_STATUS_SUCCESS) STOP "failed to initialize MKL_SPBLAS matrix"
    this%descr%type=SPARSE_MATRIX_TYPE_GENERAL 
    this%descr%diag=SPARSE_DIAG_NON_UNIT
    this%descr%mode=SPARSE_FILL_MODE_LOWER

    stat = MKL_SPARSE_CONVERT_CSR(H,SPARSE_OPERATION_NON_TRANSPOSE,this%H)
    if(stat /= 0) STOP "error setting H_ee sparse to CSR"
    stat=MKL_SPARSE_DESTROY(H)
    if(stat /= 0) STOP "error destroying H"

    Call this%init_base(lat,H_coo%op_l,H_coo%op_r)
end subroutine 


subroutine eval_single(this,E,i_m,lat)
    use m_derived_types, only: lattice
    ! input
    class(t_H_mkl_csr),intent(in)   :: this
    type(lattice), intent(in)       :: lat
    integer, intent(in)             :: i_m
    ! output
    real(8), intent(out)            :: E
    ! internal
    real(8),pointer             :: modes_l(:),modes_r(:)
    real(8),allocatable,target  :: vec_l(:),vec_r(:)
    real(C_DOUBLE)              :: tmp(this%dimH(1),1)

    integer(C_int)        :: stat
    type(SPARSE_MATRIX_T) :: vec
    real(8),external      :: ddot !blas routine

    Call lat%point_order(this%op_l,this%dimH(1),modes_l,vec_l)
    Call lat%point_order(this%op_r,this%dimH(2),modes_r,vec_r)


    !is it smarter to multiply in the other direction first with 1 entry only?-> probably easier to use different implementation 
    !where this is required
    Call create_sparse_vec(i_m,modes_r,this%dim_mode(2),this%dimH(2),vec)
    stat=mkl_sparse_d_spmmd(SPARSE_OPERATION_NON_TRANSPOSE , this%H , vec , SPARSE_LAYOUT_ROW_MAJOR , tmp ,1)
    if(stat/=SPARSE_STATUS_SUCCESS) ERROR STOP "mkl error"
    E=ddot(this%dimH(1),modes_l,1,tmp,1)

    stat=mkl_sparse_destroy(vec)
    if(stat/=SPARSE_STATUS_SUCCESS) ERROR STOP "mkl error"
    nullify(modes_l,modes_r)
    if(allocated(vec_l)) deallocate(vec_l)
    if(allocated(vec_r)) deallocate(vec_r)
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

#endif
end module