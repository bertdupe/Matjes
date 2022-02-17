module m_grp_sym
use, intrinsic :: iso_fortran_env, only : output_unit
use m_derived_types, only : t_cell
use m_vector, only : norm
use m_constants, only : pi
use m_user_info
use m_sym_public
use m_symmetry_base
use m_type_lattice, only : lattice

private
public :: find_group!,read_symmetries,get_num_sym_file,get_sym_local

contains


subroutine find_group(my_lattice,my_motif)
implicit none
type(lattice), intent(in) :: my_lattice
type(t_cell), intent(in)  :: my_motif
! internal variables
integer :: number_sym,number_sym_lat
integer, allocatable :: sym_index(:)
real(kind=8) :: time
class(pt_grp),allocatable :: base_symmetries,all_symmetries,my_symmetries
! first step determine the lattice symmetries

time=0.0d0
call user_info(output_unit,time,'calculating the symmetry operations',.false.)

call set_sym_type(base_symmetries)
call set_sym_type(my_symmetries)
call set_sym_type(all_symmetries)

call base_symmetries%init_base()
call base_symmetries%load_base()
call base_symmetries%get_all_symetries(number_sym,all_symmetries)

allocate(sym_index(number_sym),source=0)

call all_symmetries%get_latt_sym(my_lattice%areal,number_sym_lat,sym_index,my_lattice%periodic,my_lattice%dim_lat)

write(output_unit,'(/,a,I2,a,/)') 'The lattice has  ',number_sym_lat,'  symmetrie operations'

call all_symmetries%get_pt_sym(my_lattice,number_sym,sym_index,my_motif)

call my_symmetries%init(number_sym)
call my_symmetries%load(sym_index,all_symmetries)

call my_symmetries%write_sym()

deallocate(base_symmetries,all_symmetries,my_symmetries)

call user_info(output_unit,time,'done',.true.)

end subroutine find_group


end module m_grp_sym
