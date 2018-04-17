module m_lattice
! Array of spins on the lattice (matrix)
real(kind=8), allocatable,target :: Spin(:,:,:,:,:)
! Array of the table of Neighbours
Integer, allocatable :: tableNN(:,:,:,:,:,:), indexNN(:,:)
! array that contains 1 or 0 depending of the presence or absence of spins
! very important in the case of non periodic boundary conditions
Integer, allocatable :: masque(:,:,:,:)
end module m_lattice

subroutine create_lattice(mag_lattice,motif)
use m_lattice, only : spin
use m_derived_types
implicit none
type (cell), intent(in) :: motif
type (lattice), intent(inout) :: mag_lattice
! internal
integer :: dim_lat(3), nmag, Ntot
! slopes
integer :: i,j,k,l

dim_lat=mag_lattice%dim_lat
nmag=count(motif%i_mom)

Ntot=product(dim_lat)*nmag


allocate(spin(3,dim_lat(1),dim_lat(2),dim_lat(3),nmag))
spin=0.0d0
allocate(mag_lattice%l_modes(dim_lat(1),dim_lat(2),dim_lat(3),nmag))

do l=1,nmag
   do k=1,dim_lat(3)
      do j=1,dim_lat(2)
         do i=1,dim_lat(1)

       mag_lattice%l_modes(i,j,k,l)%w=>spin(:,i,j,k,l)
         enddo
      enddo
   enddo
enddo


end subroutine create_lattice
