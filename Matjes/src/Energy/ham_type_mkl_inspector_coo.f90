!submodule (m_H_type_coo) m_H_type_mkl_inspector_coo
module m_H_type_mkl_inspector_coo
#if defined(CPP_MATMUL_MKL_COO) && defined(CPP_MKL_SPBLAS)
!Hamiltonian type specifications using mkl in coordinate sparse format 
use m_H_type
use MKL_SPBLAS
use m_H_type_coo
USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_DOUBLE,C_INT
implicit none


type,extends(t_H) :: t_H_mkl_coo
    private
    integer               :: dimH(2)=0
    type(SPARSE_MATRIX_T) :: H
    TYPE(matrix_descr)    :: descr
    real(C_DOUBLE),allocatable     :: val(:)
    integer,allocatable   :: rowind(:),colind(:)
contains
    !necessary t_H routines
	procedure :: eval_single=>eval_single
	procedure :: eval_all=>eval_all
	procedure :: set_H=>set_H
    procedure :: set_H_1=> set_H_1
    procedure :: add_H=> add_H
    procedure :: destroy=> destroy
    !TOTO FINAL
end type

interface t_H_mkl_coo
	procedure :: dummy_constructor
end interface t_H_mkl_coo
 
private
public t_H,t_H_mkl_coo
contains 

subroutine destroy(this)
	class(t_H_mkl_coo),intent(inout)    :: this
    integer     ::  stat

    !inherited data
    Call this%set_prepared(.false.)
    deallocate(this%op_l,this%op_r)

    stat=mkl_sparse_destroy(this%H)
    if(stat/=SPARSE_STATUS_SUCCESS) STOP 'failed to destroy t_h_mkl_coo type in m_H_type_mkl_inspector_coo'
end subroutine

type(t_H_mkl_coo) function dummy_constructor()
	!might want some initialization for H and descr, but should work without
	!continue 
end function 

subroutine add_H(this,H_add)
	class(t_H_mkl_coo),intent(inout)    :: this
	class(t_H),intent(in)               :: H_add


    STOP "IMPLEMENT ADDIND FOR t_H_mkl_coo in m_H_type_mkl_inspector_coo if really necessary"
    !addition directly implemented in MKL_SPBLAS only for CSR/BSR format <- difficult

end subroutine 

subroutine set_H_1(this,line,Hval,Hval_ind,order,lat)
    use m_derived_types, only: lattice
	class(t_H_mkl_coo),intent(inout)    :: this

	type(lattice),intent(in)    	:: lat
    integer,intent(in)              :: order(2)
    real(8),intent(in)              :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
    integer,intent(in)              :: Hval_ind(:,:)
    integer,intent(in)              :: line(:,:)

    !local
	type(t_H_coo)           :: H_coo

    Call H_coo%set_H_1(line,Hval,Hval_ind,order,lat)
    Call set_from_Hcoo(this,H_coo)
end subroutine 

subroutine set_H(this,energy_in,lat)
    use m_derived_types, only: operator_real_order_N,lattice
	class(t_H_mkl_coo),intent(inout)    :: this
	type(operator_real_order_N)         :: energy_in
	type(lattice),intent(in)    	    :: lat

    !local
	type(t_H_coo)           :: H_coo
    integer                 :: nnz
    integer(C_int)          :: stat


    Call H_coo%set_H(energy_in,lat)
    Call set_from_Hcoo(this,H_coo)
end subroutine 

subroutine eval_single(this,E,i_m,lat)
	use m_derived_types, only: lattice
	! input
	class(t_H_mkl_coo),intent(in)    :: this
	type(lattice), intent(in) 		:: lat
	integer, intent(in) 			:: i_m
	! output
	real(kind=8), intent(out) 		:: E
    ! internal
    real(C_DOUBLE)    :: tmp(this%dimH(1))
    real(C_DOUBLE)    :: alpha,beta
    integer(C_int)    :: stat
    real(8),pointer   :: modes_l(:),modes_r(:)
    integer           :: dim_modes(2)

    Call lat%set_order_point(this%op_l(1),modes_l)
    Call lat%set_order_point(this%op_r(1),modes_r)
    dim_modes(1)=lat%get_order_dim(this%op_l(1))
    dim_modes(2)=lat%get_order_dim(this%op_r(1))

    tmp=0.0d0
    alpha=1.0d0;beta=0.0d0

    stat=mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,alpha,this%H,this%descr,modes_r,beta,tmp)
    if(stat/=SPARSE_STATUS_SUCCESS) STOP "failed MKL_SPBLAS routine in eval_single  of m_H_type_mkl_inspector_coo"
    E=dot_product(modes_l((i_m-1)*dim_modes(1)+1:i_m*dim_modes(1)),tmp((i_m-1)*dim_modes(1)+1:i_m*dim_modes(1)))
	nullify(modes_l,modes_r)
end subroutine 


subroutine eval_all(this,E,lat)
    use m_derived_types, only: lattice
	class(t_H_mkl_coo),intent(in)    :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(out)            :: E
    ! internal
    real(C_DOUBLE)    :: tmp(this%dimH(1))
    real(C_DOUBLE)    :: alpha,beta
    integer(C_int)    :: stat
    real(8),pointer   :: modes_l(:),modes_r(:)

    Call lat%set_order_point(this%op_l(1),modes_l)
    Call lat%set_order_point(this%op_r(1),modes_r)

    tmp=0.0d0
    alpha=1.0d0;beta=0.0d0
    stat=mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,alpha,this%H,this%descr,modes_r,beta,tmp)
    E=dot_product(modes_l,tmp)
    if(stat/=SPARSE_STATUS_SUCCESS) STOP "failed MKL_SPBLAS routine in eval_all of m_H_type_mkl_inspector_coo"
	nullify(modes_l,modes_r)
    
end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       PRIVATE UTILITY                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_from_Hcoo(this,H_coo)
	type(t_H_coo),intent(inout)         :: H_coo
	type(t_H_mkl_coo),intent(inout)     :: this

    !local
    integer                 :: nnz
    integer(C_int)          :: stat

    Call H_coo%pop_par(this%dimH,nnz,this%val,this%rowind,this%colind)

    stat=mkl_sparse_d_create_coo(this%H, SPARSE_INDEX_BASE_ONE , this%dimH(1) , this%dimH(2) , nnz , this%rowind , this%colind , this%val)

    if(stat/=SPARSE_STATUS_SUCCESS) STOP "failed to initialize MKL_SPBLAS matrix"
    this%descr%type=SPARSE_MATRIX_TYPE_GENERAL 
    this%descr%diag=SPARSE_DIAG_NON_UNIT
    this%descr%mode=SPARSE_FILL_MODE_LOWER

    allocate(this%op_l,source=H_coo%op_l)
    allocate(this%op_r,source=H_coo%op_r)
    Call this%set_prepared(.true.)
end subroutine




#endif
end module
