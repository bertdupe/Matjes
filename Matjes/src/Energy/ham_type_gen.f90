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
TODO 
#else
use m_H_type_manual
#endif


use m_derived_types, only : lattice
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

    function energy_all(ham,lat)result(E)
        !get all energies from an energy array
        class(t_H),intent(in)       ::  ham(:)
        class(lattice),intent(in)   ::  lat
        real(8)                     ::  E

        real(8)     ::  tmp_E(size(ham))
        integer     ::  i

        E=0.0d0
        do i=1,size(ham)
          Call ham(i)%eval_all(tmp_E(i),lat)
        enddo
        E=sum(tmp_E)
    end function
     
    function energy_single(ham,i_m,lat)result(E)
        !get all energies from an energy array
        class(t_H),intent(in)       ::  ham(:)
        integer,intent(in)          ::  i_m
        class(lattice),intent(in)   ::  lat
        real(8)                     ::  E

        real(8)     ::  tmp_E(size(ham))
        integer     ::  i

        E=0.0d0
        do i=1,size(ham)
          Call ham(i)%eval_single(tmp_E(i),i,lat)
        enddo
        E=sum(tmp_E)
    end function

end module
