       module m_measure_temp

       contains
!
! modulke that allows the temperature measurements in spin dynamics
!

! initialization of the temperature measurements
!
       subroutine init_temp_measure(check,check1,check2,check3)
       implicit none
       real(kind=8), intent(out) :: check(2),check1,check2,check3 

       check=0.0d0
       check1=0.0d0
       check2=0.0d0
       check3=0.0d0

       end subroutine init_temp_measure


!
! update the temperature
!

       subroutine update_temp_measure(check1,check2,spin,Beff)
       use m_vector, only : cross
       implicit none
       real(kind=8), intent(inout) :: check1,check2
       real(kind=8), intent(in) ::  spin(3),Beff(3)

       check1=check1+sum(cross(spin,Beff)**2)
       check2=check2+dot_product(spin,Beff)

       end subroutine update_temp_measure

!
! get the temperature
!

       subroutine get_temp(security,check,kt)
       use m_constants, only : k_B
       implicit none
       real(kind=8), intent(out) :: security(2)
       real(kind=8), intent(in) ::  kt,check(2)
 
       security=0.0d0  

       security(1)=check(1)/check(2)/2.0d0/k_B
       security(2)=(check(1)/check(2)/2.0d0-kT)/k_B*1000.0d0

       end subroutine get_temp

       end module m_measure_temp
