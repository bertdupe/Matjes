module m_H_tb_public
use m_H_tb_base
use m_H_tb_coo
use m_H_tb_dense
use m_H_tb_csr
use m_TB_types, only: parameters_TB_IO_H
public

contains
subroutine set_H(H,io)
    use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
    class(H_tb),allocatable,intent(inout)   :: H
    type(parameters_TB_IO_H),intent(in)     :: io

    if(allocated(H)) STOP "CANNOT set H which is already set"
    if(io%sparse)then
        write(output_unit,'(2/A/)') "Chose sparse feast algoritm for tight-binding Hamiltonian"
        allocate(H_feast_csr::H)
    else
        select case(io%i_diag)
        case(1)
            write(output_unit,'(2/A/)') "Chose lapack zheevd algoritm for tight-binding Hamiltonian"
            allocate(H_zheevd::H)
        case(2)
            write(output_unit,'(2/A/)') "Chose lapack zheev algoritm for tight-binding Hamiltonian"
            allocate(H_zheev::H)
        case(3)
            write(output_unit,'(2/A/)') "Chose dense feast algoritm for tight-binding Hamiltonian"
            allocate(H_feast_den::H)
        case(4)
            write(output_unit,'(2/A/)') "Chose lapack zheevr algoritm for tight-binding Hamiltonian"
            allocate(H_zheevr::H)
        case default
            write(error_unit,'(2/A,I6,A)') "Unable to choose dense tight-binding Hamiltonian as TB_diag=",io%i_diag," is not implemented"
            STOP "CHECK INPUT"
        end select
    endif
end subroutine



end module
