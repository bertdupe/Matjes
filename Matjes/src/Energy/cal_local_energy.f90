module m_local_energy
!module which sets the correct energy set/get subroutines depending on the preprocessor flags

#ifdef CPP_MATMUL_MKL_SPARSE
use m_local_energy_mkl, only: local_energy=>local_energy_mkl
use m_local_energy_mkl, only: set_E_matrix=>set_E_matrix_mkl
use m_local_energy_mkl, only: kill_E_matrix=>kill_E_matrix_mkl
use m_local_energy_mkl, only: sum_energy=>sum_energy_mkl
#elif defined CPP_MATMUL_EIGEN_SPARSE
use m_local_energy_eigen, only: local_energy=>local_energy_eigen
use m_local_energy_eigen, only: set_E_matrix=>set_E_matrix_eigen
#else
use m_local_energy_manual, only: local_energy=>local_energy_manual
use m_local_energy_manual, only: set_E_matrix=>set_E_matrix_manual
use m_local_energy_manual, only: kill_E_matrix=>kill_E_matrix_manual
use m_local_energy_manual, only: sum_energy=>sum_energy_manual
!use m_local_energy_manual, only: sum_energy=>sum_energy_manual
#endif
implicit none
public
end module m_local_energy
