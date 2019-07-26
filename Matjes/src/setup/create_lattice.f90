module m_lattice
! Array of spins on the lattice (matrix)
real(kind=8), allocatable,save,target :: modes(:,:,:,:,:)


!!!!!!
! to be deleted
! Array of the table of Neighbours
Integer, allocatable :: tableNN(:,:,:,:,:,:), indexNN(:,:)
! array that contains 1 or 0 depending of the presence or absence of spins
! very important in the case of non periodic boundary conditions
Integer, allocatable :: masque(:,:,:,:)

private
public :: modes,create_lattice,masque,tableNN
contains

subroutine create_lattice(mag_lattice,motif,ext_param)
use m_derived_types
implicit none
type(cell), intent(in) :: motif
type(lattice), intent(inout) :: mag_lattice
type(simulation_parameters), intent(in) :: ext_param
! internal
integer :: dim_lat(3), nmag, Ntot, N_dim_order_param
! slopes
integer :: i,j,k,l

dim_lat=mag_lattice%dim_lat
nmag=count(motif%atomic(:)%moment.gt.0.0d0)

N_dim_order_param=get_num_mode(motif,ext_param)


!allocate(modes(N_dim_order_param,dim_lat(1),dim_lat(2),dim_lat(3),nmag))
allocate(modes(3,dim_lat(1),dim_lat(2),dim_lat(3),nmag))
modes=0.0d0

allocate(mag_lattice%l_modes(dim_lat(1),dim_lat(2),dim_lat(3),nmag))

do l=1,nmag
   do k=1,dim_lat(3)
      do j=1,dim_lat(2)
         do i=1,dim_lat(1)

       mag_lattice%l_modes(i,j,k,l)%w=>modes(:,i,j,k,l)

         enddo
      enddo
   enddo
enddo


end subroutine create_lattice

!!!!!!!!!!!!!!!!!!!!!!!!
! function that gets the number of order parameter
!!!!!!!!!!!!!!!!!!!!!!!!
integer function get_num_mode(motif,ext_param)
use m_derived_types
implicit none
type(cell), intent(in) :: motif
type(simulation_parameters), intent(in) :: ext_param
! internal parameters
integer :: N_dim_order_param,nmag
real(kind=8) :: Field(3)

N_dim_order_param=0
nmag=count(motif%atomic(:)%moment.gt.0.0d0)
N_dim_order_param=N_dim_order_param+nmag*3
! electric field
Field=ext_param%E_ext%value
if (sqrt(Field(1)**2+Field(2)**2+Field(3)**2).gt.1.0d-8) then
   N_dim_order_param=N_dim_order_param+3
   write(6,'(a)') 'electric field found'
   write(6,'(a,2x,I5)') 'Hamiltonian has the dimension',N_dim_order_param**2
endif

! Magnetic field
Field=ext_param%H_ext%value
if (sqrt(Field(1)**2+Field(2)**2+Field(3)**2).gt.1.0d-8) then
   N_dim_order_param=N_dim_order_param+3
   write(6,'(a)') 'magnetic field found'
   write(6,'(a,2x,I5)') 'Hamiltonian has the dimension',N_dim_order_param**2
endif

! Temperature
if ((ext_param%ktini%value.gt.1.0d-8).or.(ext_param%ktfin%value.gt.1.0d-8)) then
   N_dim_order_param=N_dim_order_param+1
   write(6,'(a)') 'Temperature found'
   write(6,'(a,2x,I5)') 'Hamiltonian has the dimension',N_dim_order_param**2
endif

get_num_mode=N_dim_order_param

end function get_num_mode

end module m_lattice
