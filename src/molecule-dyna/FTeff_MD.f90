module m_FTeff_MD

use m_randist
use m_random_number_library
use m_constants, only : hbar

private
public :: draw_stocha_integrals !Bolzmann
contains


subroutine draw_stocha_integrals(sigma_u,sigma_v,c_uv,delta_ug,delta_vg)
	implicit none
   	real(8),intent(in) :: sigma_v(:),sigma_u(:), c_uv(:)  !gaussian distrib parameters
   	real(8),intent(inout) :: delta_vg(:), delta_ug(:) !stochastic integrals of u and v to draw from bivariate gaussian 
	!internal
	real(8) :: temp(2) !dummies
	integer :: i
	
	do i=1,size(delta_ug)
		!draw in normal distrib of stdev 1, mean 0
		temp(1)=randist(1.0d0)
		temp(2)=randist(1.0d0)
		write(*,*)'temp1 and 2= ',temp(1),' ',temp(2)
		!draw in bivariate gaussian distribution
		delta_ug(i) = sigma_u(i) * temp(1)
		delta_vg(i) = sigma_v(i) * (c_uv(i)*temp(1) + sqrt(1.0d0- c_uv(i)**2) *temp(2) )
		
		write(*,*)'delta_ug(i)= ',delta_ug(i),'delta_vg(i)= ',delta_vg(i)
	enddo
	
end subroutine draw_stocha_integrals







!not used
subroutine Bolzmann(kt,damping,masses,FT)
implicit none
real(kind=8), intent(in) :: kt,damping,masses(:)
real(kind=8), intent(inout) :: FT(:)
! internal
integer :: i

do i=1,size(FT)
   FT(i)=sqrt(damping*masses(i))*randist(kt)
enddo

end subroutine Bolzmann



end module m_FTeff_MD
