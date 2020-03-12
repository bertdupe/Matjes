module m_dipole_energy
use m_dipolar_field
!
! contains all the routines to calculate the dipole dipole energy interaction
!
!

private
public :: get_dipole_E

contains

!
! calculate the double sum of the dipole dipole energy
!
real(kind=8) function get_dipole_E(iomp)
use m_constants, only : mu_B
implicit none
! input
integer, intent(in) :: iomp
! internal variable
integer :: i,j,N
! internal variable
real(kind=8) :: rc(3),ss,B(3)

get_dipole_E=0.0d0
B=0.0d0

call get_dipole_B(B,iomp)

!
! the dipolar is 2 times larger to take into account the sommation on the atoms
! so the energy should be divided by 2
!

get_dipole_E=-dot_product(mode_dipole(iomp)%w,B)/2.0d0

end function get_dipole_E






!
! get the dipole dipole position matrix elements
!
!function get_distances_dipole(r1,r2,periodic)
!implicit none
!integer, intent(in) :: N
!logical, intent(in) :: periodic(3)
!real(kind=8) :: get_distances_dipole(3)
!! internal
!
!write(*,*) 'get_positions_dipole'
!
!end function get_distances_dipole

!!!! part of the FFT dipole dipole interaction
#ifndef CPP_BRUTDIP
! prepare the demag tensor
!      if (i_dip) then
!       write(6,'(a)') "preparing spatial contribution to dipole convolution "
!       call setup_dipole(dim_lat,Periodic_log,net,count(my_motif%i_mom),size(world))
!       do k=1,dim_lat(3)*count(my_motif%i_mom)
!        do j=1,dim_lat(2)
!         do i=1,dim_lat(1)
!         mmatrix(1:3,i,j,k)=spin(4:6,i,j,k,mod(k-1,count(my_motif%i_mom))+1)*spin(7,i,j,k,1)
!         enddo
!        enddo
!       enddo
!       write(6,'(a)') 'FFT dipole dipole is set up'
!       else
!       allocate(ntensor(1,1,1,1),mmatrix(1,1,1,1),ctrans(1,1,1),rtrans(1,1,1), &
!     & hcomplex(1,1,1,1),mcomplex(1,1,1,1),hreal(1,1,1,1))
!
!       mmatrix=0.0d0
!       rtrans=0.0d0
!       hreal=0.0d0
!       ctrans=dcmplx(0.d0,0.d0)
!       mcomplex=dcmplx(0.d0,0.d0)
!       hcomplex=dcmplx(0.d0,0.d0)
!       ntensor=dcmplx(0.d0,0.d0)
!      endif
!#endif
#endif

end module m_dipole_energy
