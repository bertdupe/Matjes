module m_eval_BTeff

 interface randist
    module procedure gaussianran
    module procedure wiener
 end interface

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

   write(*,*) kt
   stop

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



!!!! the different ways to generate random forces
real(kind=8) function gaussianran(kt)
implicit none
real(kind=8), intent(in) :: kt
real(kind=8) :: Choice1, Choice2

Choice1=10.0d0
Choice2=10.0d0

do while (Choice1**2+Choice2**2.gt.1.0d0)
  call RANDOM_NUMBER(Choice1)
  call RANDOM_NUMBER(Choice2)
  Choice1=(Choice1*2.0d0-1.0d0)
  Choice2=(Choice2*2.0d0-1.0d0)
enddo

!
! according to the PhD of Schieback (2010) from Constanz p. 35
!

gaussianran=sqrt(2.0d0*kT)*Choice1*sqrt(-2.0d0*dlog(Choice1**2+Choice2**2)/(Choice1**2+Choice2**2))

end function


!!!! wiener process
real(kind=8) function wiener()
use mtprng
implicit none
!dummy
type(mtprng_state) :: state
real(kind=8) :: Choice

#ifdef CPP_MRG
Choice=mtprng_rand_real1(state)
#else
CALL RANDOM_NUMBER(Choice)
#endif

wiener=(Choice*2.0d0-1.0d0)

end function wiener

end module m_eval_BTeff
