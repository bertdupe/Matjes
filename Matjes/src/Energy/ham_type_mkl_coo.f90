module  m_H_type_mkl_coo
#if defined(CPP_MATMUL_MKL_COO) && !defined(CPP_MKL_SPBLAS)
!Hamiltonian type specifications using mkl in coordinate sparse format 
use m_H_type
use m_H_type_coo


type,extends(t_H) :: t_H_mkl_coo
    private
    integer             :: dimH(2)=0 !dimension of Hamiltonian
    integer             :: nnz=0  !number of entries in sparse matrix
    real(8),allocatable :: val(:)
    integer,allocatable :: rowind(:),colind(:)
contains
    !necessary t_H routines
	procedure :: eval_single=>eval_single
	procedure :: eval_all=>eval_all
	procedure :: set_H=>set_H
	procedure :: set_H_1=>set_H_1
end type
private
public t_H,t_H_mkl_coo
contains 


subroutine set_H_1(this,line,Hval,Hval_ind,order,lat)
    use m_derived_types, only: lattice
	class(t_H_mkl_coo),intent(inout)    :: this

	type(lattice),intent(in)    	:: lat
    integer,intent(in)              :: order(2)
    real(8),intent(in)              :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
    integer,intent(in)              :: Hval_ind(:,:)
    integer,intent(in)              :: line(:,:)

	type(t_H_coo)   :: H_coo

    Call H_coo%set_H_1(line,Hval,Hval_ind,order,lat)
    Call H_coo%pop_par(this%dimH,this%nnz,this%val,this%rowind,this%colind)
    allocate(this%op_l,source=H_coo%op_L)
    allocate(this%op_r,source=H_coo%op_r)

end subroutine 




subroutine set_H(this,energy_in,lat)
    use m_derived_types, only: operator_real_order_N,lattice
	class(t_H_mkl_coo),intent(inout)    :: this
	type(operator_real_order_N)         :: energy_in
	type(lattice),intent(in)    	    :: lat

	type(t_H_coo)   :: H_coo

    Call H_coo%set_H(energy_in,lat)
    Call H_coo%pop_par(this%dimH,this%nnz,this%val,this%rowind,this%colind)
    allocate(this%op_l,source=H_coo%op_L)
    allocate(this%op_r,source=H_coo%op_r)

end subroutine 

subroutine eval_single(this,E,i_m,lat)
	use m_derived_types, only: lattice
	! input
	class(t_H_mkl_coo),intent(in)   :: this
	type(lattice), intent(in) 		:: lat
	integer, intent(in) 			:: i_m
	! output
	real(kind=8), intent(out) 		:: E
    ! internal
    real(8),allocatable             :: tmp(:)

    real(8),pointer     :: modes_l(:),modes_r(:)
    integer             :: dim_modes(2)

    if(this%dimH==0) STOP "Hamiltonian not set"
    Call lat%set_order_point(this%op_l(1),modes_l)
    Call lat%set_order_point(this%op_r(1),modes_r)
    dim_modes(1)=lat%get_order_dim(this%op_l(1))
    dim_modes(2)=lat%get_order_dim(this%op_r(1))

    allocate(tmp(dim_modes(1)))
    Call mkl_dcoogemv('N',size(modes_r),this%val,this%rowind,this%colind,this%nnz,modes_r,tmp)
    E=dot_product(modes_l((i_m-1)*dim_modes(1)+1:i_m*dim_modes(1)),tmp((i_m-1)*dim_modes(1)+1:i_m*dim_modes(1)))
end subroutine 


subroutine eval_all(this,E,lat)
    use m_derived_types, only: lattice
	class(t_H_mkl_coo),intent(in)    :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(out)            :: E
    ! internal
    real(8),allocatable             :: tmp(:)

    real(8),pointer     :: modes_l(:),modes_r(:)
    integer             :: dim_modes(2)

    if(this%dimH==0) STOP "Hamiltonian not set"
    Call lat%set_order_point(this%op_l(1),modes_l)
    Call lat%set_order_point(this%op_r(1),modes_r)
    dim_modes(1)=lat%get_order_dim(this%op_l(1))
    dim_modes(2)=lat%get_order_dim(this%op_r(1))
    allocate(tmp(size(modes_r)))
    Call mkl_dcoogemv('N',size(modes_r),this%val,this%rowind,this%colind,this%nnz,modes_r,tmp)
    E=dot_product(modes_l,tmp)
    nullify(modes_l,modes_r)
    deallocate(tmp)
    
end subroutine 

#endif
end module
