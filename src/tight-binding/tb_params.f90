module m_tb_params
use m_tb_types
implicit none
public

type(parameters_TB),public,protected :: TB_params

contains
subroutine set_TB_params(Ncell)
    use m_rw_TB, only: get_parameters_io_TB
    integer,intent(in)      ::  Ncell

    Call get_parameters_io_TB(TB_params)

    if(TB_params%io_H%nb_spin == 1)then
        TB_params%is_mag=.False.
    elseif (TB_params%io_H%nb_spin == 2)then
        TB_params%is_mag=.True.
    else 
        STOP "Unexpected value for TB_params%io_H%nb_spin"
    endif

    if(any(TB_params%io_H%delta/=0.0d0)) TB_params%is_sc=.True.

    if(TB_params%is_sc .and. .not. TB_params%is_mag) STOP "Trying to use SC without using the magnetic part"

    !set dimensions of the Hamiltonian
    if(TB_params%is_sc) TB_params%H%nsc=2
    if(TB_params%is_mag) TB_params%H%nspin=2
    TB_params%H%ncell=Ncell
    TB_params%H%norb=TB_params%io_H%nb_orbitals
    Call TB_params%H%upd()

    TB_params%H%sparse=TB_params%io_H%sparse
    TB_params%H%i_diag=TB_params%io_H%i_diag
    TB_params%H%estNe=TB_params%io_H%estNe
    TB_params%H%extE=TB_params%io_H%extE
    TB_params%H%diag_acc=TB_params%io_H%diag_acc
    if(TB_params%H%extE(1)>=TB_params%H%extE(2))then
        TB_params%H%extE(1)=-9.0d99!Huge(1.0d0)
        TB_params%H%extE(2)=9.0d99!Huge(1.0d0)
    endif

    !TB_params%H%rearrange=TB_params%io_H%rearrange !this is not really implementend and most probably breakes dos or something else
end subroutine



end module
