module m_H_public
!module used to choose the correct version of treating the Hamiltonian 
!and provide all necessary types for other routines

#if defined CPP_MKL_SPBLAS
use m_H_sparse_mkl
#elif defined CPP_EIGEN_H
use m_H_eigen
use m_H_eigen_mem
#elif defined CPP_DENSE && defined(CPP_BLAS)
use m_H_dense_blas
#elif defined CPP_DENSE
use m_H_dense
#else
!TODO MORE
use m_H_manual
#endif

use m_derived_types, only : lattice
implicit none
public

contains
    subroutine get_Htype(H_out)
        class(t_h),intent(out),allocatable      :: H_out 
#if defined CPP_MKL_SPBLAS
       allocate(t_H_mkl_csr::H_out)
#elif defined CPP_EIGEN_H
       !allocate(t_H_eigen::H_out)
       allocate(t_H_eigen_mem::H_out)   !mem saves transpose H as well which accelerates some operations, but costs memeory (add option)
#elif defined CPP_DENSE && defined(CPP_BLAS)
       allocate(t_H_dense_blas::H_out)
#elif defined CPP_DENSE
       allocate(t_H_dense::H_out)
#else
       allocate(t_H_manual::H_out)
#endif
    end subroutine
    
    subroutine get_Htype_N(H_out,N)
        class(t_h),intent(out),allocatable      :: H_out (:)
        integer,intent(in)                      :: N
#if defined CPP_MKL_SPBLAS
       allocate(t_H_mkl_csr::H_out(N))
#elif defined CPP_EIGEN_H
       !allocate(t_H_eigen::H_out(N))
       allocate(t_H_eigen_mem::H_out(N))    !mem saves transpose H as well which accelerates some operations, but costs memeory (add option)
#elif defined CPP_DENSE && defined(CPP_BLAS)
       allocate(t_H_dense_blas::H_out(N))
#elif defined CPP_DENSE
       allocate(t_H_dense::H_out(N))
#else
       allocate(t_H_manual::H_out(N))
#endif
    end subroutine

    subroutine Bcast_Harr(ham,comm)
        use mpi_basic                
        class(t_H),intent(inout),allocatable    ::  ham(:)
        type(mpi_type),intent(in)               ::  comm
#ifdef CPP_MPI
        integer     ::  i,N

        if(comm%ismas) N=size(ham)
        Call MPI_Bcast(N,1, MPI_INTEGER, comm%mas, comm%com,i)
        if(.not.comm%ismas) Call get_Htype_N(ham,N)
        do i=1,size(ham)
            Call Ham(i)%bcast(comm)
        enddo
#else
        continue
#endif
    end subroutine

    function energy_all(ham,lat)result(E)
        !get all energies from an energy array
        class(t_H),intent(in)       ::  ham(:)
        class(lattice),intent(in)   ::  lat
        real(8)                     ::  E

        real(8)     ::  tmp_E(size(ham))
        
        Call energy_resolved(ham,lat,tmp_E)
        E=sum(tmp_E)
    end function

    subroutine energy_resolved(ham,lat,E)
        !get contribution-resolved energies from an energy array
        class(t_H),intent(in)       ::  ham(:)
        class(lattice),intent(in)   ::  lat
        real(8),intent(out)         ::  E(size(ham))

        integer     ::  i

        E=0.0d0
        do i=1,size(ham)
            Call ham(i)%eval_all(E(i),lat)
        enddo
    end subroutine
     
    function energy_single(ham,i_m,dim_bnd,lat)result(E)
        use m_derived_types, only: number_different_order_parameters
        !get all energies from an energy array
        class(t_H),intent(in)       :: ham(:)
        integer,intent(in)          :: i_m
        class(lattice),intent(in)   :: lat
        integer, intent(in)         :: dim_bnd(2,number_different_order_parameters)  
        real(8)                     :: E

        real(8)     ::  tmp_E(size(ham))
        integer     ::  i

        E=0.0d0
        do i=1,size(ham)
            Call ham(i)%eval_single(tmp_E(i),i_m,dim_bnd,lat)
        enddo
        tmp_E=tmp_E*real(ham(:)%mult_M_single,8)
        E=sum(tmp_E)
    end function

end module
