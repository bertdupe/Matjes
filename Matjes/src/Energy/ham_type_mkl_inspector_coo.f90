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
    integer               :: dimH=0
    type(SPARSE_MATRIX_T) :: H
    TYPE(matrix_descr)    :: descr
    real(C_DOUBLE),allocatable     :: val(:)
    integer,allocatable     :: rowind(:),colind(:)
contains
    !necessary t_H routines
	procedure :: eval_single=>eval_single
	procedure :: eval_all=>eval_all
	procedure :: set_H=>set_H
end type

interface t_H_mkl_coo
	procedure :: dummy_constructor
end interface t_H_mkl_coo
 
private
public t_H,t_H_mkl_coo
contains 

type(t_H_mkl_coo) function dummy_constructor()
	!might want some initialization for H and descr, but should work without
	!continue 
end function 


subroutine set_H(this,energy_in,lat)
    use m_derived_types, only: operator_real_order_N,lattice
	class(t_H_mkl_coo),intent(inout)    :: this
	type(operator_real_order_N)         :: energy_in
	type(lattice),intent(in)    	    :: lat

    !local
	type(t_H_coo)           :: H_coo
    integer                 :: dimH
    integer                 :: nnz
    integer(C_int)          :: stat


    Call H_coo%set_H(energy_in,lat)
    Call H_coo%pop_par(dimH,nnz,this%val,this%rowind,this%colind)
    this%dimH=dimH
    stat=mkl_sparse_d_create_coo(this%H, SPARSE_INDEX_BASE_ONE , dimH , dimH , nnz , this%rowind , this%colind , this%val)

    if(stat/=SPARSE_STATUS_SUCCESS) STOP "failed to initialize MKL_SPBLAS matrix"
    this%descr%type=SPARSE_MATRIX_TYPE_GENERAL 
    this%descr%diag=SPARSE_DIAG_NON_UNIT
    this%descr%mode=SPARSE_FILL_MODE_LOWER

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
    real(C_DOUBLE)      :: tmp(this%dimH)
    real(C_DOUBLE)      :: alpha,beta
    integer(C_int)      :: stat

    tmp=0.0d0
    alpha=1.0d0;beta=0.0d0

    stat=mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,alpha,this%H,this%descr,lat%ordpar%all_modes,beta,tmp)
    if(stat/=SPARSE_STATUS_SUCCESS) STOP "failed MKL_SPBLAS routine in eval_single  of m_H_type_mkl_inspector_coo"
    E=dot_product(lat%ordpar%all_modes((i_m-1)*lat%dim_mode+1:i_m*lat%dim_mode),tmp((i_m-1)*lat%dim_mode+1:i_m*lat%dim_mode))
end subroutine 


subroutine eval_all(this,E,lat)
    use m_derived_types, only: lattice
	class(t_H_mkl_coo),intent(in)    :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(out)            :: E
    ! internal
    real(C_DOUBLE)      :: tmp(this%dimH)
    real(C_DOUBLE)      :: alpha,beta
    integer(C_int)      :: stat

    tmp=0.0d0
    alpha=1.0d0;beta=0.0d0
    stat=mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,alpha,this%H,this%descr,lat%ordpar%all_modes,beta,tmp)
    E=dot_product(lat%ordpar%all_modes,tmp)
    if(stat/=SPARSE_STATUS_SUCCESS) STOP "failed MKL_SPBLAS routine in eval_all  of m_H_type_mkl_inspector_coo"
    
end subroutine 

#endif
end module
