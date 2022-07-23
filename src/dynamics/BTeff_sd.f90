module m_eval_BTeff
use m_random_base
use m_constants, only : hbar
use m_vector, only : norm,cross
implicit none

private
public :: langevin_bath,wiener_bath




contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! random force centered on 0 with a gaussian distribution
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine langevin_bath(is_set,damping,mode,DT,random_numbers)
logical, intent(in) :: is_set
real(kind=8), intent(in) :: mode(:,:),damping,random_numbers(:,:)
real(kind=8), intent(inout) :: DT(:,:)
! internal
real(kind=8) :: step_T(3),BT(3,size(mode,2))
integer :: i,Nspin

if (.not.is_set) return
if ((size(mode,2).ne.size(random_numbers,2)).and.(size(mode,2).ne.size(DT,2))) stop 'error in langevin_bath: size do not match'

BT=sqrt(damping*hbar)*random_numbers

Nspin=size(mode,2)

do i=1,Nspin
   step_T=cross(mode(:,i),BT(:,i),1,3)
   DT(:,i)=(step_T+damping*cross(mode(:,i),step_T,1,3))/(1+damping**2)
enddo

end subroutine langevin_bath

subroutine wiener_bath(is_set,damping,mode,DT,random_numbers)
logical, intent(in) :: is_set
real(kind=8), intent(in) :: damping,mode(:,:),random_numbers(:,:)
real(kind=8), intent(inout) :: DT(:,:)
! internal
real(kind=8) :: step_T(3),BT(3,size(mode,2))
integer :: i,Nspin

if (.not.is_set) return
if ((size(mode,2).ne.size(random_numbers,2)).and.(size(mode,2).ne.size(DT,2))) stop 'error in langevin_bath: size do not match'

Nspin=size(mode,2)

BT=random_numbers

do i=1,i
   step_T=cross(mode(:,i),BT(:,i),1,3)
   DT(:,i)=step_T+damping*cross(mode(:,i),step_T,1,3)
enddo
end subroutine wiener_bath

end module m_eval_BTeff
