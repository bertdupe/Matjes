module m_sym_public
use m_symmetry_base
use m_sym_dense
use m_sym_symlib_ccp4
!use mpi_util ,only : mpi_type, bcast
implicit none

private
public :: set_sym_type

contains

subroutine set_sym_type(my_pt_grp)
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
   class(pt_grp),allocatable   :: my_pt_grp
!   type(mpi_type),intent(in)   :: comm

   allocate(pt_grp_dense::my_pt_grp)
   write(output_unit,'(//A/)') "Using the internal symmetry routines"

end subroutine

end module
