set(main_f90 ${main_f90}
src/Hamiltonian_type/ham_public.f90
src/Hamiltonian_type/ham_base.f90
src/Hamiltonian_type/ham_coo.f90
src/Hamiltonian_type/ham_coo_based.f90
src/Hamiltonian_type/ham_cusparse.f90
src/Hamiltonian_type/ham_dense.f90
src/Hamiltonian_type/ham_dense_blas.f90
src/Hamiltonian_type/ham_deriv.f90
src/Hamiltonian_type/ham_eigen.f90
src/Hamiltonian_type/ham_eigen_mem.f90
src/Hamiltonian_type/ham_mkl_csr.f90
src/Hamiltonian_type/ham_mkl_csr_mem.f90
src/Hamiltonian_type/hamiltonian_collection.f90
src/Hamiltonian_type/fft_ham_fftw.f90
src/Hamiltonian_type/fft_ham_cufft.f90
src/Hamiltonian_type/fft_ham_base.f90
src/Hamiltonian_type/fft_ham_public.f90
src/Hamiltonian_type/fft_ham_internal.f90
src/Hamiltonian_type/work_ham_single.f90
)

set(eigen_H ${eigen_H}
src/Hamiltonian_type/eigen_H_interface.f90
src/Hamiltonian_type/eigen_H.cpp
)

set(cuda_fft ${cuda_fft}
src/Hamiltonian_type/cuda_fft_interface.f90
src/Hamiltonian_type/cuda_fft.cpp
src/Hamiltonian_type/cuda_fft.cu
)

set(cuda_H ${cuda_H}
src/Hamiltonian_type/cuda_H.cu
src/Hamiltonian_type/cuda_H_interface.f90
src/Hamiltonian_type/cuda_H_vec.cu
src/Hamiltonian_type/cuda_H_vec_interface.f90
)
