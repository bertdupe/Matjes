module m_bcast_global
!broadcasts the globally defined variables and changeable module variables
use mpi_basic
use mpi_util

contains

subroutine bcast_global_var(com)
    !collective routine which broadcasts different module and global parameters
    use m_fft_H_public, only: bcast_fft_H_mode
    use m_H_public, only: bcast_H_mode
    class(mpi_type),intent(in)  :: com
#ifdef CPP_MPI    
    Call bcast_H_mode(com)
    Call bcast_fft_H_mode(com)
#else
    continue
#endif
end subroutine




end module

