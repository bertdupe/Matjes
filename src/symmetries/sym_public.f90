module m_sym_public
use m_symmetry_base
use m_sym_dense
use m_sym_symlib_ccp4
use m_sym_dense_smart
use m_spglib
!use mpi_util ,only : mpi_type, bcast
implicit none

integer :: mode
real(8) :: tolerance

private
public :: set_sym_type,set_mode

contains

subroutine set_sym_type(my_pt_grp)
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
   class(pt_grp),allocatable   :: my_pt_grp
!   type(mpi_type),intent(in)   :: comm


   select case (mode)
     case(1)
       allocate(pt_grp_dense::my_pt_grp)
       write(output_unit,'(//A/)') "Using the internal symmetry routines; method Bertrand"

     case(2)
       allocate(pt_grp_dense_smart::my_pt_grp)
       write(output_unit,'(//A/)') "Using the internal symmetry routines; method Melanie (smart method)"

     case(3)
#ifdef CPP_SPGLIB
       allocate(sym_spglib::my_pt_grp)
       write(output_unit,'(//A/)') "Using the spglib symmetry libraries"
#else
       write(output_unit,'(//A/)') "spglib symmetry libraries will not be used"
#endif

     case(4)
  !     allocate(sym_spglib::my_pt_grp)
  !     write(output_unit,'(//A/)') "Using the spglib symmetry libraries"
  stop "not coded yet"

     case default
       stop 'Symmetry mode not implemented'
   end select

   my_pt_grp%tol_sym=tolerance

end subroutine

subroutine set_mode(mode_in,tolerance_in)
integer, intent(in) :: mode_in
real(8), intent(in) :: tolerance_in


mode=mode_in
tolerance=tolerance_in

end subroutine

end module
