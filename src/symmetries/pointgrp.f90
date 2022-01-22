module m_grp_sym
use, intrinsic :: iso_fortran_env, only : output_unit
use m_derived_types, only : t_cell
use m_vector, only : norm
use m_constants, only : pi
use m_user_info
use m_sym_public
use m_symmetry_base

private
public :: find_group!,read_symmetries,get_num_sym_file,get_sym_local

contains


subroutine find_group(areal,my_motif,periodic,dim_lat)
implicit none
real(kind=8), intent(in) :: areal(3,3)
type(t_cell), intent(in) :: my_motif
logical, intent(in) :: periodic(:)
integer, intent(in) :: dim_lat(:)
! internal variables
integer :: number_sym,number_sym_lat
integer, allocatable :: sym_index(:)
real(kind=8) :: time
class(pt_grp),allocatable :: all_symmetries,my_symmetries
! first step determine the lattice symmetries

time=0.0d0
call user_info(output_unit,time,'calculating the symmetry operations',.false.)

call set_sym_type(all_symmetries)
call set_sym_type(my_symmetries)

call all_symmetries%init_base()
call all_symmetries%load_base()

number_sym_lat=all_symmetries%get_N_sym()
number_sym=number_sym_lat
allocate(sym_index(number_sym_lat),source=0)

call all_symmetries%get_latt_sym(areal,number_sym_lat,sym_index)

write(output_unit,'(/,a,I2,a,/)') 'The lattice has  ',number_sym_lat,'  symmetrie operations'

call all_symmetries%get_pt_sym(areal,number_sym,sym_index,my_motif,periodic,dim_lat)

call my_symmetries%init(number_sym)
call my_symmetries%load(sym_index,all_symmetries)

call my_symmetries%write_sym()

deallocate(all_symmetries,my_symmetries)

call user_info(output_unit,time,'done',.true.)

end subroutine find_group



end module m_grp_sym
