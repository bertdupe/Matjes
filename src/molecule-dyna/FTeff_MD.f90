module m_FTeff_MD

use m_randist
use m_random_number_library
use m_constants, only : hbar
use m_random_number_library
use m_random_public

private
public :: draw_stocha_integrals 
contains


subroutine draw_stocha_integrals(thermal_noise,rand1,rand2,sigma_u,sigma_v,c_uv,delta_ug,delta_vg)
	implicit none
	class(ranbase), intent(inout) :: thermal_noise !for drawing N*dim(mode) random numbers
	real(8),intent(inout) :: rand1(:),rand2(:) !to draw random numbers
   real(8),intent(in) :: sigma_v(:),sigma_u(:), c_uv(:)  !gaussian distrib parameters
   real(8),intent(inout) :: delta_vg(:), delta_ug(:) !stochastic integrals of u and v to draw from bivariate gaussian 
	
	!this could be replaced by a single call to VSL_RNG_METHOD_GAUSSIANMV_BOXMULLER2
	call thermal_noise%get_extract_list(0.0d0,1.0d0,rand1) !mean,sigma,resu
	call thermal_noise%get_extract_list(0.0d0,1.0d0,rand2) !mean,sigma,resu

	delta_ug = sigma_u *rand1
	delta_vg = sigma_v * ( c_uv*rand1 + sqrt(1.0d0 - c_uv**2) * rand2 )

end subroutine draw_stocha_integrals



end module m_FTeff_MD
