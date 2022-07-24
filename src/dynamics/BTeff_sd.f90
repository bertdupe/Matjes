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
subroutine langevin_bath(is_set,kt,damping,mode,DT,random_numbers,time_step)
logical, intent(in) :: is_set
real(kind=8), intent(in) :: mode(:,:),damping,random_numbers(:,:),kt,time_step
real(kind=8), intent(inout) :: DT(:,:)
! internal
real(kind=8) :: step_T(3),BT(3,size(mode,2))
integer :: i,Nspin

Nspin=size(mode,2)

if (.not.is_set) return
if ((Nspin.ne.size(random_numbers,2)).or.(Nspin.ne.size(DT,2))) stop 'error in langevin_bath: size do not match'

BT=sqrt(2.0d0*damping*hbar*kT/time_step)*random_numbers

do i=1,Nspin
   step_T=cross(mode(:,i),BT(:,i))
   DT(:,i)=-(step_T+damping*cross(mode(:,i),step_T))/(1.0d0+damping*damping)
enddo

end subroutine langevin_bath

subroutine wiener_bath(is_set,kt,damping,mode,DT,random_numbers,time_step)
logical, intent(in) :: is_set
real(kind=8), intent(in) :: damping,mode(:,:),random_numbers(:,:),kt,time_step
real(kind=8), intent(inout) :: DT(:,:)
! internal
real(kind=8) :: step_T(3),BT(3,size(mode,2))
integer :: i,Nspin

if (.not.is_set) return
if ((size(mode,2).ne.size(random_numbers,2)).and.(size(mode,2).ne.size(DT,2))) stop 'error in langevin_bath: size do not match'

Nspin=size(mode,2)

BT=random_numbers

do i=1,i
   step_T=cross(mode(:,i),BT(:,i))
   DT(:,i)=step_T+damping*cross(mode(:,i),step_T)
enddo
end subroutine wiener_bath

end module m_eval_BTeff
