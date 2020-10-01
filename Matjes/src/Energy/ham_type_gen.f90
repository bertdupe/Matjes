module m_Htype_gen
!module used to choose the correct version of treating the Hamiltonian 
!and provide all necessary types for other routines

#if defined   CPP_MATMUL_MKL_CSR && defined(CPP_MKL_SPBLAS)
use m_H_type_mkl_inspector_csr
#elif defined CPP_MATMUL_MKL_COO && defined(CPP_MKL_SPBLAS)
use m_H_type_mkl_inspector_coo
#elif defined CPP_MATMUL_MKL_CSR 
use m_H_type_mkl_csr
#elif defined CPP_MATMUL_MKL_COO
use m_H_type_mkl_coo
#elif defined CPP_MATMUL_EIGEN_SPARSE
TOTO 
#else
use m_H_type_manual
#endif
implicit none
public

contains
subroutine get_Htype(H_out)
    class(t_h),intent(out),allocatable      :: H_out 
#if defined CPP_MATMUL_MKL_CSR
   allocate(H_out,source=t_H_mkl_csr())
#elif defined CPP_MATMUL_MKL_COO
   allocate(H_out,source=t_H_mkl_coo())
!#elif defined CPP_MATMUL_EIGEN_SPARSE
!    TOTO 
#else
   allocate(H_out,source=t_H_manual())
#endif
end subroutine

subroutine get_Htype_N(H_out,N)
    class(t_h),intent(out),allocatable      :: H_out (:)
    integer,intent(in)                      :: N
#if defined CPP_MATMUL_MKL_CSR
   allocate(H_out(N),source=t_H_mkl_csr())
#elif defined CPP_MATMUL_MKL_COO
   allocate(H_out(N),source=t_H_mkl_coo())
!#elif defined CPP_MATMUL_EIGEN_SPARSE
!    TOTO 
#else
   allocate(H_out(N),source=t_H_manual())
#endif
end subroutine


end module
