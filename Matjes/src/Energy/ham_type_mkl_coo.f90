module m_H_type_mkl_coo
#if defined CPP_MATMUL_MKL_COO
!Hamiltonian type specifications using mkl in coordinate sparse format 
use m_H_type
use m_H_type_coo


type,extends(t_H) :: t_H_mkl_coo
    private
    integer             :: dimH=0 !dimension of Hamiltonian
    integer             :: nnz=0  !number of entries in sparse matrix
    real(8),allocatable :: val(:)
    integer,allocatable :: rowind(:),colind(:)
contains
    !necessary t_H routines
	procedure :: eval_single=>eval_single
	procedure :: eval_all=>eval_all
	procedure :: set_H=>set_H
end type
private
public t_H,t_H_mkl_coo
contains 

subroutine set_H(this,energy_in,lat)
    use m_derived_types, only: operator_real_order_N,lattice
	class(t_H_mkl_coo),intent(inout)    :: this
	type(operator_real_order_N)         :: energy_in
	type(lattice),intent(in)    	    :: lat

	type(t_H_coo)   :: H_coo

    Call H_coo%set_H(energy_in,lat)
    Call H_coo%pop_par(this%dimH,this%nnz,this%val,this%rowind,this%colind)

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
    real(8)             :: tmp(this%dimH)

    if(this%dimH==0) STOP "Hamiltonian not set"

    Call mkl_dcoogemv('N',this%dimH,this%val,this%rowind,this%colind,this%nnz,lat%ordpar%all_modes,tmp)
    E=dot_product(lat%ordpar%all_modes((i_m-1)*lat%dim_mode+1:i_m*lat%dim_mode),tmp((i_m-1)*lat%dim_mode+1:i_m*lat%dim_mode))
end subroutine 


subroutine eval_all(this,E,lat)
    use m_derived_types, only: lattice
	class(t_H_mkl_coo),intent(in)    :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(out)            :: E
    ! internal
    real(8)             :: tmp(this%dimH)

    if(this%dimH==0) STOP "Hamiltonian not set"
    
    Call mkl_dcoogemv('N',this%dimH,this%val,this%rowind,this%colind,this%nnz,lat%ordpar%all_modes,tmp)
    E=dot_product(lat%ordpar%all_modes,tmp)
    
end subroutine 

#endif
end module
