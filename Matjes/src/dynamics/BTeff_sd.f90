module m_eval_BTeff
use m_randist

private
public :: langevin_bath,wiener_bath
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! random force centered on 0 with a gaussian distribution
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine langevin_bath(kt,damping,mode,BT,DT,size_mode)
use m_constants, only : hbar
use m_vector, only : norm,cross
implicit none
integer, intent(in) :: size_mode
real(kind=8), intent(in) :: kt,mode(:),damping
real(kind=8), intent(inout) :: BT(:),DT(:)
! internal
real(kind=8) :: step_T(3)
integer :: i,j

do i=1,size_mode
   BT(i)=sqrt(damping*hbar)*randist(kt)
enddo

do i=1,size_mode/3
   j=3*(i-1)+1
   step_T=cross(mode(j:j+2),BT(j:j+2),1,3)
   DT(j:j+2)=(step_T(:)+damping*cross(mode(j:j+2),step_T(j:j+2),1,3))/(1+damping**2)
enddo

end subroutine langevin_bath

subroutine wiener_bath(kt,damping,mode,BT,DT,size_mode)
use m_vector, only : norm,cross
implicit none
!logical, intent(in) :: stmtemp
integer, intent(in) :: size_mode
real(kind=8), intent(in) :: kt,damping,mode(:)
real(kind=8), intent(inout) :: BT(:),DT(:)
! internal
real(kind=8) :: step_T(3)
integer :: i,j

do i=1,size_mode
   BT(i)=randist()
enddo

do i=1,size_mode/3
   j=3*(i-1)+1
   step_T=cross(mode(j:j+2),BT(j:j+2),1,3)
   DT(j:j+2)=step_T(:)+damping*cross(mode(j:j+2),step_T(j:j+2),1,3)
enddo
end subroutine wiener_bath

end module m_eval_BTeff
