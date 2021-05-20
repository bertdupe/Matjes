module m_H_public
!module used to choose the correct version of treating the Hamiltonian 
!and provide all necessary types for other routines

use m_H_sparse_mkl
use m_H_cusparse
use m_H_eigen
use m_H_eigen_mem
use m_H_dense_blas
use m_H_dense

use m_derived_types, only : lattice
implicit none

public
integer,private ::  mode=-1 !stores which implementation is used (1:mkl, 2: eigen_mem, 3: eigen, 4:dense_blas, 5: dense, 6: cuda, -1: not initialized)

contains
    subroutine set_Ham_mode_io(fname_in)
        use m_io_files_utils, only : open_file_read,close_file
        use m_io_utils,only: get_parameter
        character(*),intent(in),optional    :: fname_in
        !internal
        character(*),parameter              :: fname_default='input'
        character(:), allocatable           :: fname

        integer ::  io_param
#if defined CPP_CUDA
        mode=6
#elif defined CPP_MKL
        mode=1
#elif defined CPP_EIGEN
        mode=2
#elif defined CPP_BLAS
        mode=4
#else
        mode=5
#endif

        if(present(fname_in))then
            fname=fname_in
        else
            fname=fname_default
        endif
        io_param=open_file_read(fname)
        call get_parameter(io_param,fname,'Hamiltonian_mode',mode)
        call close_file(fname,io_param)

        select case(mode)
        case(1)
            write(*,'(/A/)') "Using Hamiltonian implementation: MKL Sparse"
        case(2)
            write(*,'(/A/)') "Using Hamiltonian implementation: Eigen (saving transpose as well)"
        case(3)
            write(*,'(/A/)') "Using Hamiltonian implementation: Eigen"
        case(4)
            write(*,'(/A/)') "Using Hamiltonian implementation: Dense using blas"
        case(5)
            write(*,'(/A/)') "Using Hamiltonian implementation: Dense"
        case(6)
            write(*,'(/A/)') "Using Hamiltonian implementation: Cuda"
        case default
            ERROR STOP "READ IN UNIMPLEMENTED Hamiltonian_mode"
        end select
    end subroutine

    subroutine get_Htype(H_out)
        class(t_h),intent(out),allocatable      :: H_out 
        select case(mode)
        case(1)
#if defined CPP_MKL
            allocate(t_H_mkl_csr::H_out)
#else
            ERROR STOP "CANNOT USE t_H_mkl_csr Sparse Hamiltonian  implementation without sparse mkl (CPP_MKL)"
#endif
        case(2)
#if defined CPP_EIGEN
            allocate(t_H_eigen_mem::H_out)
#else
            ERROR STOP "CANNOT USE t_H_eigen_mem Sparse Hamiltonian  implementation without Eigen (CPP_EIGEN)"
#endif
        case(3)
#if defined CPP_EIGEN
            allocate(t_H_eigen::H_out)
#else
            ERROR STOP "CANNOT USE t_H_eigen Sparse Hamiltonian implementation without Eigen (CPP_EIGEN)"
#endif
        case(4)
#if defined CPP_BLAS
            allocate(t_H_dense_blas::H_out)
#else
            ERROR STOP "CANNOT USE t_H_dense_blas Sparse Hamiltonian implementation without BLAS (CPP_BLAS)"
#endif
        case(5)
            allocate(t_H_dense::H_out)
        case(6)
#if defined CPP_CUDA
            allocate(t_H_cusparse::H_out)
#else
            ERROR STOP "CANNOT USE t_H_cusparse Sparse Hamiltonian implementation without Cuda (CPP_CUDA)"
#endif
        case(-1)
            ERROR STOP "Cannot allocate Hamiltonian type, mode has not been set using set_H_mode_io"
        case default
            ERROR STOP "UNEXPECTED MODE SET DECIDING HAMILTONIAN IMPLEMENTATION"
        end select
    end subroutine
    
    subroutine get_Htype_N(H_out,N)
        class(t_h),intent(out),allocatable      :: H_out (:)
        integer,intent(in)                      :: N

        select case(mode)
        case(1)
#if defined CPP_MKL
            allocate(t_H_mkl_csr::H_out(N))
#else
            ERROR STOP "CANNOT USE t_H_mkl_csr Sparse Hamiltonian  implementation without sparse mkl (CPP_MKL)"
#endif
        case(2)
#if defined CPP_EIGEN
            allocate(t_H_eigen_mem::H_out(N))
#else
            ERROR STOP "CANNOT USE t_H_eigen_mem Sparse Hamiltonian  implementation without Eigen (CPP_EIGEN)"
#endif
        case(3)
#if defined CPP_EIGEN
            allocate(t_H_eigen::H_out(N))
#else
            ERROR STOP "CANNOT USE t_H_eigen Sparse Hamiltonian implementation without Eigen (CPP_EIGEN)"
#endif
        case(4)
#if defined CPP_BLAS
            allocate(t_H_dense_blas::H_out(N))
#else
            ERROR STOP "CANNOT USE t_H_dense_blas Hamiltonian implementation without BLAS (CPP_BLAS)"
#endif
        case(5)
            allocate(t_H_dense::H_out(N))
        case(6)
#if defined CPP_CUDA
            allocate(t_H_cusparse::H_out(N))
#else
            ERROR STOP "CANNOT USE t_H_cusparse Sparse Hamiltonian implementation without Cuda (CPP_CUDA)"
#endif
        case(-1)
            ERROR STOP "Cannot allocate Hamiltonian type, mode has not been set using set_H_mode_io"
        case default
            ERROR STOP "UNEXPECTED MODE SET DECIDING HAMILTONIAN IMPLEMENTATION"
        end select
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

    subroutine bcast_H_mode(com)
        use mpi_util
        class(mpi_type),intent(in)  :: com

        Call bcast(mode,com)
    end subroutine

end module
