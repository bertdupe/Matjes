module m_measure_temp
use m_constants, only : k_B
use m_vector, only : cross,norm_cross
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit

implicit none

private
public  :: get_temp_measure

contains
!
! module that allows the temperature measurements in spin dynamics
!


!
! update the temperature
!

subroutine get_temp_measure(spin,Beff,T_measured)
real(kind=8), intent(in)  :: spin(:,:),Beff(:,:)
real(kind=8), intent(out) :: T_measured

integer :: i,Nspin
real(8) :: check1,check2

Nspin=size(Beff,2)
check1=0.0d0
check2=0.0d0

if (Nspin.ne.size(spin,2)) STOP "Beff and spin do not have the same size"

do i=1,Nspin
   check1=check1+norm_cross(spin(:,i),Beff(:,i))**2
   check2=check2+dot_product(spin(:,i),Beff(:,i))
enddo

T_measured=check1/check2/2.0d0/k_B

write(output_unit,'(a,2x,f16.6)') 'Temperature (K)', T_measured

end subroutine get_temp_measure

end module m_measure_temp
