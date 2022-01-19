module m_FT_Ham_public
!use m_FT_Ham_coo
!use m_FT_csr
use m_parameters_FT_Ham
use m_H_type
use m_FT_ham_base
use m_FT_Ham_zheev
use m_FT_Ham_dense
implicit none

private :: set_H_single
public :: set_H

interface set_H
    module procedure set_H_single
end interface

contains

subroutine set_H_single(H,io)
    use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
    class(FT_Ham_base),allocatable,intent(inout)      :: H
    type(parameters_FT_HAM_IO), intent(in)            :: io
    logical,save    :: nsaid=.true.

    if(allocated(H)) STOP "CANNOT set H which is already set"
!      select case(mode)
!      case(1)
!#ifdef CPP_MKL
!        write(output_unit,'(2/A/)') "Chose sparse feast algorithm for FT Hamiltonian"
!        allocate(FT_H_feast_csr::H)
!#else
!        write(error_unit,'(//A)') "CANNOT use sparse feast algorithm without CPP_MKL"
!        write(error_unit,'(A)')   "Disable sparse FT desription or compile with USE_MKL"
!        ERROR STOP
!#endif
!      end select
!    else
        select case(io%i_diag)
        case(1)
#ifdef CPP_LAPACK
            if(nsaid) write(output_unit,'(2/A/)') "Choose lapack zheevd algoritm for FT Hamiltonian"
            allocate(FT_Ham_zheev::H)
#else
            write(error_unit,'(//A)') "CANNOT use lapack diagonalization algorithm without CPP_LAPACK"
            ERROR STOP
#endif
        case(2)
            if(nsaid) write(output_unit,'(2/A/)') "Choose internal diagonalization algoritm for FT Hamiltonian"
            allocate(FT_H_dense::H)
!        case(3)
!#ifdef CPP_MKL
!            if(nsaid) write(output_unit,'(2/A/)') "Chose dense feast algoritm for tight-binding Hamiltonian"
!            allocate(H_feast_den::H)
!#else
!            write(error_unit,'(//A)') "CANNOT use dense feast algorithm without CPP_MKL"
!            write(error_unit,'(A)')   "Choose different TB_diag  or compile with USE_MKL"
!            ERROR STOP
!#endif
!        case(4)
!#ifdef CPP_LAPACK
!            if(nsaid) write(output_unit,'(2/A/)') "Chose lapack zheevr algoritm for tight-binding Hamiltonian"
!            allocate(H_zheevr::H)
!#else
!            write(error_unit,'(//A)') "CANNOT use lapack diagonalization algorithm without CPP_LAPACK"
!            ERROR STOP
!#endif
!        case default
!            write(error_unit,'(2/A,I6,A)') "Unable to choose dense tight-binding Hamiltonian as TB_diag=",io%i_diag," is not implemented"
!            STOP "CHECK INPUT"
        end select
!    endif
    nsaid=.false.
end subroutine

end module m_FT_Ham_public
