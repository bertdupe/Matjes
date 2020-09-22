#ifdef CPP_MATMUL_MKL_SPARSE

module m_local_energy_mkl
use m_basic_types, only : vec_point
use m_derived_types, only : point_shell_Operator,lattice
use m_modes_variables, only : point_shell_mode
use m_local_energy_manual, only: get_matrix_sparse_manual
implicit none

#ifdef CPPMKL_CSR
real(8),allocatable :: val_csr(:)
integer,allocatable :: i_csr(:),j_csr(:)
#else
real(8),allocatable :: val(:)
integer,allocatable :: rowind(:),colind(:)
#endif
integer             :: dimH
!public :: set_H_sparse,energy_sparse,dimH


private
public :: local_energy_mkl,set_E_matrix_mkl,kill_E_matrix_mkl
public :: sum_energy_mkl

contains

subroutine set_E_matrix_mkl(dim_mode)
    integer,intent(in)      ::  dim_mode
    integer                 ::  nnz
    integer                 ::  tmp,info
    integer,parameter       ::  job(8)=[2,1,1,0,0,0,0,0]
#ifdef CPPMKL_CSR
    real(8),allocatable :: val(:)
    integer,allocatable :: rowind(:),colind(:)
#endif
    external mkl_dcsrcoo

    Call get_matrix_sparse_manual(dim_mode,val,rowind,colind,dimH)
#ifdef CPPMKL_CSR
    nnz=size(val)
    allocate(val_csr(nnz),source=0.0d0)
    allocate(j_csr(nnz),source=0)
    allocate(i_csr(dimH+1),source=0)
    tmp=0
    info=0
    Call mkl_dcsrcoo(job,dimH,val_csr,j_csr,i_csr,nnz,val,rowind,colind,tmp,info)
    deallocate(val,rowind,colind)
#endif
end subroutine


subroutine local_energy_mkl(E,iomp,lat)
    !evaluates energy for a single cell spin modification
    real(8),intent(out) :: E
    type(lattice),intent(in)  :: lat
    integer,intent(in)  :: iomp

    integer             :: dimH
    integer             :: nnz
    integer             :: N_site
    real(8),allocatable :: tmp(:)
    integer             :: dim_mode
    external mkl_dcoogemv
  
    dim_mode=lat%dim_mode
    dimH=product(lat%dim_lat)*dim_mode
    allocate(tmp(dimH))
#ifdef CPPMKL_CSR
    nnz=size(val_csr)
    Call mkl_dcsrgemv('N',dimH,val_csr,i_csr,j_csr,lat%ordpar%all_modes,tmp)
#else
    nnz=size(val)
    Call mkl_dcoogemv('N',dimH,val,rowind,colind,nnz,lat%ordpar%all_modes,tmp)
#endif
    E=dot_product(lat%ordpar%all_modes((iomp-1)*dim_mode+1:iomp*dim_mode),tmp((iomp-1)*dim_mode+1:iomp*dim_mode))
end subroutine


subroutine sum_energy_mkl(E,lat)
    real(8),intent(out)         :: E
    type(lattice),intent(in)    :: lat
    !internal
    integer             :: nnz
    integer             :: N_site
    real(8)             :: tmp(dimH)
    external mkl_dcoogemv
   
#ifdef CPPMKL_CSR
    nnz=size(val_csr)
    Call mkl_dcsrgemv('N',dimH,val_csr,i_csr,j_csr,lat%ordpar%all_modes,tmp)
#else
    nnz=size(val)
    Call mkl_dcoogemv('N',dimH,val,rowind,colind,nnz,lat%ordpar%all_modes,tmp)
#endif
    E=dot_product(tmp,lat%ordpar%all_modes)
end subroutine

subroutine kill_E_matrix_mkl()

#ifdef CPPMKL_CSR
    deallocate(val_csr,i_csr,j_csr)
#else
    deallocate(val,rowind,colind)
#endif
    write(6,'(a)') 'Energy matrix deallocated'
end subroutine 

end module
#endif
