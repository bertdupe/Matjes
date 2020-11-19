module m_projection

public

contains

! Convert effective fields to the format used in VPO
! project a vector v on a vector u and reconstruct the vector in the base of the u's
subroutine project_force(v,u,fxyz)
implicit none
real(kind=8), intent(in) :: v(:),u(:)
real(kind=8), intent(out) :: fxyz(:)
! internal variables
real(kind=8) :: tmp

! coordinate of v on u
tmp = dot_product(v,u)

! there might be a factor 2 missing here (to check)
fxyz = v - tmp*u

end subroutine project_force

end module m_projection
