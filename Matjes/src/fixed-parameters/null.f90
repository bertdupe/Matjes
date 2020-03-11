module m_null
use m_derived_types, only : site_Ham
!
! create a matrix of 0 to points to
!
type(site_Ham), target, public, protected, save :: nunull

private
public :: get_null_matrix

contains

subroutine get_null_matrix(dim_ham)
implicit none
integer, intent(in) :: dim_ham

allocate(nunull%H(dim_ham,dim_ham))

nunull%H=0.0d0

end subroutine get_null_matrix

end module m_null
