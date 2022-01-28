module m_FTeff_MD

use m_randist
use m_random_number_library
use m_constants, only : hbar
use m_random_number_library

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
	!write(*,*)'size(delta_ug)',size(delta_ug),' shape(delta_ug)=',shape(delta_ug)
	
	!write to file
	!open(1,file='normal.dat')
	!open(2,file='delta_ug_vg.dat')
	do i=1,size(delta_ug)
		!draw in normal distrib of stdev 1, mean 0
		!temp(1)=randist(1.0d0)
		!temp(2)=randist(1.0d0)

		temp(1)=rand_normal(0.0d0,1.0d0)
		temp(2)=rand_normal(0.0d0,1.0d0)
	!write(1,*) temp(:)

		!draw in bivariate gaussian distribution
	!	u = m1 + s1 * X;
	!	v = m2 + s2 * (rho * X + sqrt(1 - rho^2) * Y);

		delta_ug(i) = sigma_u(i) * temp(1)
		delta_vg(i) = sigma_v(i) * ( c_uv(i)*temp(1) + sqrt(1.0d0 - c_uv(i)**2) * temp(2) )
	 !	write(2,*) delta_ug(i), '	',delta_vg(i)
		!write(*,*)'delta_ug(i)= ',delta_ug(i),'delta_vg(i)= ',delta_vg(i)
	enddo
	!close(1)
	!close(2)
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
