module m_sym_dense
use m_symmetry_base
implicit none

private
public :: pt_grp_dense

type, extends(pt_grp) ::pt_grp_dense
contains
 procedure :: init

end type

contains

subroutine init(this,N)
 type(pt_grp), intent(out) :: this
 integer     , intent(in)  :: N

 integer :: i

 allocate(this(N))
 ! set to zero
 do i=1,N
     this(i1)%mat = 0.d0
     this(i1)%name=""
 end do
 this%n_sym=N

end subroutine

end module
