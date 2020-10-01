module m_H_type_mkl_csr
#if defined(CPP_MATMUL_MKL_CSR) && !defined(CPP_MKL_SPBLAS)
!Hamiltonian type specifications using mkl in coordinate sparse format 
use m_H_type
use m_H_type_coo, only: t_h_coo

type,extends(t_H) :: t_H_mkl_csr
    private
    integer             :: dimH=0 !dimension of Hamiltonian
    integer             :: nnz=0  !number of entries in sparse matrix
    real(8),allocatable :: val_csr(:)
    integer,allocatable :: i_csr(:),j_csr(:)
contains
	procedure :: eval_single=>eval_single
	procedure :: eval_all=>eval_all
	procedure :: set_H=>set_H
end type
private
public t_H,t_H_mkl_csr
contains 

subroutine set_H(this,energy_in,lat)
    use m_derived_types, only: operator_real_order_N,lattice
	!need to get rid of dim_mode input
	class(t_H_mkl_csr),intent(inout)  :: this
	type(operator_real_order_N)       :: energy_in
	type(lattice),intent(in)          :: lat

	type(t_H_coo)   :: H_coo

    real(8),allocatable :: val(:)
    integer,allocatable :: rowind(:),colind(:)
    integer             :: tmp,info
    integer,parameter   :: job(8)=[2,1,1,0,0,0,0,0]
    external mkl_dcsrcoo

    Call H_coo%set_H(energy_in,lat) 
    Call H_coo%pop_par(this%dimH,this%nnz,val,rowind,colind)

    allocate(this%val_csr(this%nnz),source=0.0d0)
    allocate(this%j_csr(this%nnz),source=0)
    allocate(this%i_csr(this%dimH+1),source=0)
    tmp=0
    info=0
    Call mkl_dcsrcoo(job,this%dimH,this%val_csr,this%j_csr,this%i_csr,this%nnz,val,rowind,colind,tmp,info)
    deallocate(val,rowind,colind)


end subroutine 

subroutine eval_single(this,E,i_m,lat)
	use m_derived_types, only: lattice
	! input
	class(t_H_mkl_csr),intent(in)   :: this
	type(lattice), intent(in) 		:: lat
	integer, intent(in) 			:: i_m
	! output
	real(kind=8), intent(out) 		:: E
    ! internal
    real(8)       :: tmp(this%dimH)

    Call mkl_dcsrgemv('N',this%dimH,this%val_csr,this%i_csr,this%j_csr,lat%ordpar%all_modes,tmp)
    E=dot_product(lat%ordpar%all_modes((i_m-1)*lat%dim_mode+1:i_m*lat%dim_mode),tmp((i_m-1)*lat%dim_mode+1:i_m*lat%dim_mode))

end subroutine 


subroutine eval_all(this,E,lat)
    use m_derived_types, only: lattice
	class(t_H_mkl_csr),intent(in)   :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(out)            :: E
    ! internal
    real(8)       :: tmp(this%dimH)
    
    Call mkl_dcsrgemv('N',this%dimH,this%val_csr,this%i_csr,this%j_csr,lat%ordpar%all_modes,tmp)
    E=dot_product(lat%ordpar%all_modes,tmp)
    
end subroutine 

#endif
end module
