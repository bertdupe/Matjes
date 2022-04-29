module m_3spin
use, intrinsic :: iso_fortran_env, only : output_unit
use m_input_H_types, only: io_H_sp3
use m_derived_types, only: lattice
use m_neighbor_pair, only: pair_dat_t, get_pair_dat_M

implicit none
character(len=*),parameter  :: ham_desc="3-spin interaction"
private
public :: read_sp3_input, get_3spin

contains

subroutine read_sp3_input(io_param,fname,io)
    use, intrinsic :: iso_fortran_env, only : output_unit,error_unit
    integer,intent(in)              :: io_param
    character(len=*), intent(in)    :: fname
    type(io_H_sp3),intent(out)      :: io
    character(len=*),parameter      :: var_name="M_3spin"


end subroutine

subroutine get_3spin(Ham,io,lat)
    !get 4-spin interaction in t_H Hamiltonian format
    use m_H_public, only: t_H, get_Htype
    use m_mode_public
    use m_coo_mat

    class(t_H),intent(inout)    :: Ham
    type(io_H_sp3),intent(in)   :: io
    type(lattice),intent(in)    :: lat

end subroutine

end module
