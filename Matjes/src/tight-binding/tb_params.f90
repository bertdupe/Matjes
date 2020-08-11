module m_tb_params
use m_tb_types
implicit none
public

type(parameters_TB),public,protected :: TB_params

contains
subroutine set_TB_params(Ncell,TB_pos_ext)
    use m_rw_TB, only: get_parameters_io_TB
    integer,intent(in)      ::  Ncell
    integer,intent(in)      ::  TB_pos_ext(2)

    Call get_parameters_io_TB(TB_params)

    if(TB_params%io_H%nb_spin == 1)then
        TB_params%is_mag=.False.
    elseif (TB_params%io_H%nb_spin == 2)then
        TB_params%is_mag=.True.
    else 
        STOP "Unexpected value for TB_params%io_H%nb_spin"
    endif

    if(any(TB_params%io_H%delta/=0.0d0))then
#ifdef CPP_SC 
        TB_params%is_sc=.True.
#else
        STOP "trying to use superconductivity but CPP_SC is not set" 
#endif
    endif

    if(TB_params%is_sc .and. .not. TB_params%is_mag) STOP "Trying to use SC without using the magnetic part"

    !set dimensions of the Hamiltonian
    if(TB_params%is_sc) TB_params%H%nsc=2
    if(TB_params%is_mag) TB_params%H%nspin=2
    TB_params%H%ncell=Ncell
    TB_params%H%norb=TB_params%io_H%nb_orbitals
    TB_params%H%pos_ext=TB_pos_ext
    Call TB_params%H%upd()
end subroutine



end module
