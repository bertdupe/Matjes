       module m_lattice
! Array of spins on the lattice (matrix)
       real(kind=8), allocatable :: Spin(:,:,:,:,:)
! Array of the table of Neighbours
       Integer, allocatable :: tableNN(:,:,:,:,:,:), indexNN(:,:)
! array that contains 1 or 0 depending of the presence or absence of spins
! very important in the case of non periodic boundary conditions
       Integer, allocatable :: masque(:,:,:,:)
       end module m_lattice

       subroutine create_lattice(dim_lat,motif)
       use m_lattice, only : spin
       use m_derived_types
       implicit none
       integer, intent(in) :: dim_lat(3)
       type (cell), intent(in) :: motif
       integer :: i

       allocate(spin(7,dim_lat(1),dim_lat(2),dim_lat(3),count(motif%i_m)))
       spin=0.0d0

       end subroutine create_lattice
