module m_eval_BTeff
use m_random_base
use m_constants, only : hbar
use m_vector, only : norm,cross,cross_NM,normalize
implicit none

private
public :: langevin_bath,wiener_bath


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! random force centered on 0 with a gaussian distribution
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine langevin_bath(is_set,kt,damping,mode,DT,random_numbers)
logical, intent(in) :: is_set
real(kind=8), intent(in) :: mode(:,:),damping,random_numbers(:,:),kt
real(kind=8), intent(out) :: DT(:,:)
! internal
real(kind=8) :: step_FL(3,size(mode,2)),BT(3,size(mode,2)),step_DL(3,size(mode,2))
real(8)     ::  M_norm(size(mode,1),size(mode,2))
integer :: i,Nspin

   step_FL=0.0d0
   step_DL=0.0d0
   DT=0.0d0

   Nspin=size(mode,2)
   if (.not.is_set) return
   if ((Nspin.ne.size(random_numbers,2)).or.(Nspin.ne.size(DT,2))) stop 'error in langevin_bath: size do not match'

   BT=random_numbers

   M_norm=mode
   Call normalize(M_norm)
   call cross_NM(M_norm,BT,step_FL)
   step_DL=BT+damping*step_FL
   call cross_NM(M_norm,step_DL,DT)

   DT=-DT/(1.0d0+damping*damping)

end subroutine langevin_bath

subroutine wiener_bath(is_set,kt,damping,mode,DT,random_numbers)
logical, intent(in) :: is_set
real(kind=8), intent(in) :: damping,mode(:,:),random_numbers(:,:),kt
real(kind=8), intent(out) :: DT(:,:)
! internal
real(kind=8) :: step_T(3,size(mode,2)),BT(3,size(mode,2))
integer :: i,Nspin

DT=0.0d0

if (.not.is_set) return
if ((size(mode,2).ne.size(random_numbers,2)).and.(size(mode,2).ne.size(DT,2))) stop 'error in langevin_bath: size do not match'

Nspin=size(mode,2)

   BT=random_numbers
   call cross_NM(mode,BT,DT)
   call cross_NM(mode,DT,step_T)

   DT=-(DT+damping*step_T)

end subroutine wiener_bath

end module m_eval_BTeff
