module m_eval_BTeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module that calculate the different random fields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

interface calculate_BTeff
  module procedure random_field
end interface

private
public :: calculate_BTeff
contains

subroutine random_field(stmtemp,kt,BT)
use m_randist
implicit none
logical, intent(in) :: stmtemp
real(kind=8), intent(in) :: kt
real(kind=8), intent(inout) :: BT(:)
! internal
real(kind=8) :: BT_norm

BT_norm=sqrt(BT(1)**2+BT(2)**2+BT(3)**2)

if (kt.gt.1.0d-10) then
  if (stmtemp) then
!    steptemp=cross(spini,(/randist(kt),randist(kt),randist(kt)/))*htor(i_x,i_y,i_z)/maxh
  else
    if (BT_norm.lt.1.0d-8) BT=(/randist(kt),randist(kt),randist(kt)/)
  endif
endif

end subroutine random_field

end module m_eval_BTeff
