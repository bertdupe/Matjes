module m_tb_params
use m_tb_types
implicit none
public

type(parameters_TB),public,protected :: TB_params

contains
subroutine set_TB_params()
    use m_rw_TB, only: get_parameters_io_TB
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

end subroutine



end module
