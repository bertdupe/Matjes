module m_tb_k_public
use m_tb_k_base
use m_tb_k_dense
use m_tb_k_zheev
use m_tb_k_zheevr
use m_tb_k_zheevd
use m_tb_k_zheevx
use m_tb_k_feast
use m_work_ham_single, only: work_ham
use mpi_util ,only : mpi_type, bcast
private set_H_single
interface set_Hk
    module procedure set_H_single
end interface

contains
subroutine set_H_single(Hk,mode_in,comm)
    use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
    class(H_k_base),allocatable :: Hk
    integer,intent(in)          :: mode_in
    type(mpi_type),intent(in)   :: comm
    integer     ::  mode

    mode=mode_in
    Call bcast(mode,comm)
    select case (mode)
    case (1)
#ifdef CPP_LAPACK
        allocate(H_k_zheevd::Hk)
        if(comm%ismas) write(output_unit,'(//A/)') "Set Hk-type to: zheevd"
#else
        if(comm%ismas)then
            write(error_unit,'(//A)') "ERROR, chose to use zheevd implementation of H_k_base, but lapack is not set"
            write(error_unit,'(A)') "Use different mode without lapack or compile with lapack."
        endif
        Call comm%barrier()
        ERROR STOP "CHECK input"
#endif
    case (2)
#ifdef CPP_LAPACK
        allocate(H_k_zheev::Hk)
        if(comm%ismas) write(output_unit,'(//A/)') "Set Hk-type to: zheev"
#else
        if(comm%ismas)then
            write(error_unit,'(//A)') "ERROR, chose to use zheev implementation of H_k_base, but lapack is not set"
            write(error_unit,'(A)') "Use different mode without lapack or compile with lapack."
        endif
        Call comm%barrier()
        ERROR STOP "CHECK input"
#endif
    case (3)
#ifdef CPP_MKL
        allocate(H_k_feast::Hk)
        if(comm%ismas) write(output_unit,'(//A/)') "Set Hk-type to: feast_dense"
#else
        if(comm%ismas)then
            write(error_unit,'(//A)') "ERROR, chose to use feast_dense implementation of H_k_base, but lapack is not set"
            write(error_unit,'(A)') "Use different mode without lapack or compile with lapack."
        endif
        Call comm%barrier()
        ERROR STOP "CHECK input"
#endif
    case (4)
#ifdef CPP_LAPACK
        allocate(H_k_zheevr::Hk)
        if(comm%ismas) write(output_unit,'(//A/)') "Set Hk-type to: zheevr"
#else
        if(comm%ismas)then
            write(error_unit,'(//A)') "ERROR, chose to use zheevr implementation of H_k_base, but lapack is not set"
            write(error_unit,'(A)') "Use different mode without lapack or compile with lapack."
        endif
        Call comm%barrier()
        ERROR STOP "CHECK input"
#endif

    case (5)
#ifdef CPP_LAPACK
        allocate(H_k_zheevx::Hk)
        if(comm%ismas) write(output_unit,'(//A/)') "Set Hk-type to: zheevx"
#else
        if(comm%ismas)then
            write(error_unit,'(//A)') "ERROR, chose to use zheevx implementation of H_k_base, but lapack is not set"
            write(error_unit,'(A)') "Use different mode without lapack or compile with lapack."
        endif
        Call comm%barrier()
        ERROR STOP "CHECK input"
#endif
    case default
        if(comm%ismas)then
            write(error_unit,'(//A)') "Failed to allocate H_k_base, unimlemented mode selected"
            write(error_unit,'(A,I12)') "  mode=",mode
        endif
        Call comm%barrier()
        ERROR STOP "CHECK input"
    end select
end subroutine

    

end module
