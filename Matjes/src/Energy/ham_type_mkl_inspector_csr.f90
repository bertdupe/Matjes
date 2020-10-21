
module m_H_type_mkl_inspector_csr
#if defined(CPP_MATMUL_MKL_CSR) && defined(CPP_MKL_SPBLAS)
!Hamiltonian type specifications using MKL_SPARSE inspector mkl in csr 
use m_H_type
use MKL_SPBLAS
use m_H_type_coo
use m_derived_types, only: lattice
USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_DOUBLE,C_INT


type,extends(t_H) :: t_H_mkl_csr
    private
    integer               :: dimH(2)=0
    type(SPARSE_MATRIX_T) :: H
    TYPE(matrix_descr)    :: descr
contains
    !necessary t_H routines
    procedure :: eval_single
    procedure :: eval_all   
    procedure :: set_H      
    procedure :: set_H_1    
    procedure :: set_H_mult_2   
    procedure :: destroy    
    procedure :: add_H      
    procedure :: copy       
    procedure :: optimize
    procedure :: mult_r,mult_l
    procedure :: mult_r_red,mult_l_red
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


subroutine mult_r_red(this,lat,res,op_keep)
    !multiply out right side and reduce to only keep operator corresponding to op_keep
    use m_derived_types, only: lattice
    class(t_H_mkl_csr),intent(in)   :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)   !result matrix-vector product
    integer,intent(in)              :: op_keep
    ! internal
    real(8),allocatable             :: tmp(:)   !multipied, but not reduced

    allocate(tmp(this%dimH(1)))
    Call this%mult_r(lat,tmp)
    Call lat%reduce(tmp,this%op_l,op_keep,res)
    deallocate(tmp)
end subroutine 

subroutine mult_l_red(this,lat,res,op_keep)
    !multiply out right side and reduce to only keep operator corresponding to op_keep
    use m_derived_types, only: lattice
    class(t_H_mkl_csr),intent(in)   :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)   !result matrix-vector product
    integer,intent(in)              :: op_keep
    ! internal
    real(8),allocatable             :: tmp(:)   !multipied, but not reduced

    allocate(tmp(this%dimH(2)))
    Call this%mult_l(lat,tmp)
    Call lat%reduce(tmp,this%op_r,op_keep,res)
    deallocate(tmp)
end subroutine 




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

    Call init_mode_onsite(this%op_r,this%dimH(2),lat,modes,vec)

    if(size(res)/=this%dimH(1)) STOP "size of vec is wrong"
    res=0.0d0
    stat=mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,alpha,this%H,this%descr,modes,beta,res)
    if(stat/=SPARSE_STATUS_SUCCESS) STOP "failed MKL_SPBLAS routine in mult_r of m_H_type_mkl_inspector_csr"
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

    Call init_mode_onsite(this%op_l,this%dimH(1),lat,modes,vec)

    if(size(res)/=this%dimH(2)) STOP "size of vec is wrong"
    res=0.0d0
    stat=mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE,alpha,this%H,this%descr,modes,beta,res)
    if(stat/=SPARSE_STATUS_SUCCESS) STOP "failed MKL_SPBLAS routine in mult_l of m_H_type_mkl_inspector_csr"
    
end subroutine 


subroutine optimize(this)
    class(t_H_mkl_csr),intent(inout)   :: this

    integer     :: stat

    stat=mkl_sparse_set_dotmv_hint ( this%H , SPARSE_OPERATION_NON_TRANSPOSE , this%descr , 100000 )
    stat=mkl_sparse_optimize(this%H)
    
end subroutine

subroutine set_H_1(this,line,Hval,Hval_ind,order,lat)
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
    Call H_coo%set_H_1(line,Hval,Hval_ind,order,lat)
    Call set_from_Hcoo(this,H_coo)
end subroutine 


subroutine set_H_mult_2(this,connect,Hval,Hval_ind,op_l,op_r,lat)
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
    Call H_coo%set_H_mult_2(connect,Hval,Hval_ind,op_l,op_r,lat)
    Call set_from_Hcoo(this,H_coo)
end subroutine 


subroutine copy(this,Hout)
    class(t_H_mkl_csr),intent(in)   :: this
    class(t_H),intent(inout)        :: Hout
    integer         ::  stat
    
    select type(Hout)
    class is(t_H_mkl_csr)
        Call Hout%destroy()
        Call this%copy_base(Hout)
        Hout%dimH=this%dimH
        stat=mkl_sparse_copy(this%H,this%descr,Hout%H )
        Hout%descr=this%descr
        if(stat/=SPARSE_STATUS_SUCCESS) STOP 'failed to copy Sparse_Matrix_T in m_H_type_mkl_inspector_csr'
    class default
        STOP "Cannot copy t_h_mkl_csr type with Hamiltonian that is not a class of t_h_mkl_csr"
    end select
end subroutine

subroutine add_H(this,H_add)
    class(t_H_mkl_csr),intent(inout)    :: this
    class(t_H),intent(in)               :: H_add
    
    type(SPARSE_MATRIX_T)       :: tmp_H
    integer                     :: stat
    real(C_DOUBLE),parameter    :: alpha=1.0d0


    select type(H_add)
    class is(t_H_mkl_csr)
        if(this%is_set())then
            if(.not.all(this%op_l==H_add%op_l)) STOP "CANNOT ADD hamiltonians with different op_l"
            if(.not.all(this%op_r==H_add%op_r)) STOP "CANNOT ADD hamiltonians with different op_r"
            tmp_H=this%H
            stat=mkl_sparse_d_add(sparse_operation_non_transpose,tmp_H,alpha,h_add%h,this%H)
            if(stat/=SPARSE_STATUS_SUCCESS) STOP "add_H failed in mkl_inspector_csr"
            stat=mkl_sparse_destroy(tmp_H)
            if(stat/=SPARSE_STATUS_SUCCESS) STOP 'failed to destroy t_h_mkl_coo type in m_H_type_mkl_inspector_csr'
        else
            Call H_add%copy(this)
        endif
    class default
        STOP "Cannot add t_h_mkl_csr type with Hamiltonian that is not a class of t_h_mkl_csr"
    end select

end subroutine 

subroutine destroy(this)
    class(t_H_mkl_csr),intent(inout)    :: this
    integer     ::  stat

    if(this%is_set())then
        Call this%destroy_base()
        stat=mkl_sparse_destroy(this%H)
        if(stat/=SPARSE_STATUS_SUCCESS) STOP 'failed to destroy t_h_mkl_csr type in m_H_type_mkl_inspector_csr'
        this%dimH=0
    endif
end subroutine

subroutine set_H(this,energy_in,lat)
    use m_derived_types, only: operator_real_order_N,lattice
    class(t_H_mkl_csr),intent(inout)    :: this
    type(operator_real_order_N)         :: energy_in
    type(lattice),intent(in)            :: lat

    !local
    type(t_H_coo)           :: H_coo

    if(this%is_set()) STOP "cannot set hamiltonian as it is already set"
    Call H_coo%set_H(energy_in,lat)
    Call set_from_Hcoo(this,H_coo)

end subroutine 

subroutine set_from_Hcoo(this,H_coo)
    type(t_H_coo),intent(inout)         :: H_coo
    type(t_H_mkl_csr),intent(inout)     :: this

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

    allocate(this%op_l,source=H_coo%op_l)
    allocate(this%op_r,source=H_coo%op_r)
    Call this%set_prepared(.true.)
end subroutine 


subroutine eval_single(this,E,i_m,lat)
    use m_derived_types, only: lattice
    ! input
    class(t_H_mkl_csr),intent(in)    :: this
    type(lattice), intent(in)       :: lat
    integer, intent(in)             :: i_m
    ! output
    real(kind=8), intent(out)       :: E
    ! internal
    real(C_DOUBLE)      :: tmp(this%dimH(1))
    real(8),pointer             :: modes_l(:)
    real(8),allocatable,target  :: vec_l(:)
    integer             :: dim_modes(2)

    Call init_mode_onsite(this%op_l,this%dimH(1),lat,modes_l,vec_l)
    Call this%mult_r(lat,tmp)

    STOP "UPDATE EVAL_SINGLE FOR HIGHER RANKS, AND ALSO DIM_MODES SEEMS BROKEN..." !there is some function to get those already
    ! try sparse matrix product to substitute sparse matrix times sparse vector product
    E=dot_product(modes_l((i_m-1)*dim_modes(1)+1:i_m*dim_modes(1)),tmp((i_m-1)*dim_modes(1)+1:i_m*dim_modes(1)))
end subroutine 


subroutine eval_all(this,E,lat)
    class(t_H_mkl_csr),intent(in)   :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(out)            :: E
    ! internal
    real(C_DOUBLE)             :: tmp(this%dimH(1))
    real(8),pointer            :: modes_l(:)
    real(8),allocatable,target :: vec_l(:)

    Call init_mode_onsite(this%op_l,this%dimH(1),lat,modes_l,vec_l)

    Call this%mult_r(lat,tmp)
    E=dot_product(modes_l,tmp)

    if(allocated(vec_l)) deallocate(vec_l) 
end subroutine 

subroutine init_mode_onsite(op,dimH,lat,modes,vec)
    !Subroutine that points modes to the order parameter vector according to op and dimH input
    !If size(op)>1 (i.e. dimension is folded from higher rank) allocates vec, sets it correctly
    !, and points modes
    !This only works if the unfolded order paramters are only considered on the same site
    integer,intent(in)                      :: op(:)
    integer,intent(in)                      :: dimH
    type(lattice), intent(in)               :: lat
    real(8),pointer,intent(out)             :: modes(:)
    real(8),allocatable,target,intent(out)  :: vec(:)

    if(size(op)==1)then
        Call lat%set_order_point(op(1),modes)
    else
        allocate(vec(dimH),source=0.0d0)
        Call lat%set_order_comb(op,vec)
        modes=>vec
    endif
end subroutine

#endif
end module
