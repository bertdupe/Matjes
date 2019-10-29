module m_lattice
use m_basic_types, only : var_name
use m_derived_types, only : order_parameter

! array that stores what are the order parameter used
type(order_parameter),protected,allocatable,save,target :: my_order_parameters(:)
! Array of spins on the lattice (matrix)
real(kind=8),allocatable,save,target :: modes(:,:,:,:,:)


!!!!!!
! to be deleted
! Array of the table of Neighbours
!Integer, allocatable :: tableNN(:,:,:,:,:,:), indexNN(:,:)
! array that contains 1 or 0 depending of the presence or absence of spins
! very important in the case of non periodic boundary conditions
!Integer, allocatable :: masque(:,:,:,:)

private
public :: modes,create_lattice,my_order_parameters
contains

subroutine create_lattice(mag_lattice,motif,ext_param)
use m_derived_types
implicit none
type(cell), intent(in) :: motif
type(lattice), intent(inout) :: mag_lattice
type(simulation_parameters), intent(in) :: ext_param
! internal
integer :: dim_lat(3), nmag, N_dim_order_param,N_mode
! slopes
integer :: i,j,k,l

dim_lat=mag_lattice%dim_lat
!nmag=count(motif%atomic(:)%moment.gt.0.0d0)
nmag=1

N_mode=get_num_mode(motif,ext_param)
allocate(my_order_parameters(N_mode))
! initialize the data
do i=1,N_mode
    my_order_parameters(i)%name=''
    my_order_parameters(i)%start=-1
    my_order_parameters(i)%end=-1
enddo

N_dim_order_param=get_dim_mode(motif,ext_param)

!!!!!!!!!!!!!!!!!!!!!!
!
! Create a structure that stores the order parameter of arbitraty dimanesions
! 3 if there is only the magnetization
! 6 of there is a magnetic field
! ......
!
!!!!!!!!!!!!!!!!!!!!!!
allocate(modes(N_dim_order_param,dim_lat(1),dim_lat(2),dim_lat(3),nmag))

!allocate(modes(3,dim_lat(1),dim_lat(2),dim_lat(3),nmag))
modes=0.0d0
mag_lattice%dim_mode=N_dim_order_param

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
integer function get_dim_mode(motif,ext_param)
use m_derived_types
implicit none
type(cell), intent(in) :: motif
type(simulation_parameters), intent(in) :: ext_param
! internal parameters
integer :: N_dim_order_param,nmag,i,N_mode
real(kind=8) :: Field(3)

N_dim_order_param=0
nmag=count(motif%atomic(:)%moment.gt.0.0d0)
N_dim_order_param=N_dim_order_param+nmag*3

! for each order parameter, organize the matrix of the mode
N_mode=0
do i=1,nmag
   N_mode=N_mode+1
   call order_mode(N_mode,'magnetic')
enddo

! electric field
Field=ext_param%E_ext%value
if (sqrt(Field(1)**2+Field(2)**2+Field(3)**2).gt.1.0d-8) then
   N_dim_order_param=N_dim_order_param+3
   N_mode=N_mode+1
   call order_mode(N_mode,'Efield')
   write(6,'(a)') 'electric field found'
   write(6,'(a,2x,I5)') 'Hamiltonian has the dimension',N_dim_order_param**2
endif

! Magnetic field
Field=ext_param%H_ext%value
if (sqrt(Field(1)**2+Field(2)**2+Field(3)**2).gt.1.0d-8) then
   N_dim_order_param=N_dim_order_param+3
   N_mode=N_mode+1
   call order_mode(N_mode,'Bfield')
   write(6,'(a)') 'magnetic field found'
   write(6,'(a,2x,I5)') 'Hamiltonian has the dimension',N_dim_order_param**2
endif

! Temperature
if ((ext_param%ktini%value.gt.1.0d-8).or.(ext_param%ktfin%value.gt.1.0d-8)) then
   N_dim_order_param=N_dim_order_param+1
   N_mode=N_mode+1
   call order_mode(N_mode,'temperature')
   write(6,'(a)') 'Temperature found'
   write(6,'(a,2x,I5)') 'Hamiltonian has the dimension',N_dim_order_param**2
endif

get_dim_mode=N_dim_order_param

end function get_dim_mode

!!!!!!!!!!!!!!!!!!!!!!!!
! function that gets the number of different order parameter
!!!!!!!!!!!!!!!!!!!!!!!!
integer function get_num_mode(motif,ext_param)
use m_derived_types
implicit none
type(cell), intent(in) :: motif
type(simulation_parameters), intent(in) :: ext_param
! internal parameters
integer :: N_mode,nmag
real(kind=8) :: Field(3)

N_mode=0
nmag=count(motif%atomic(:)%moment.gt.0.0d0)
N_mode=N_mode+nmag
! electric field
Field=ext_param%E_ext%value
if (sqrt(Field(1)**2+Field(2)**2+Field(3)**2).gt.1.0d-8) then
   N_mode=N_mode+1
endif

! Magnetic field
Field=ext_param%H_ext%value
if (sqrt(Field(1)**2+Field(2)**2+Field(3)**2).gt.1.0d-8) then
   N_mode=N_mode+1
endif

! Temperature
if ((ext_param%ktini%value.gt.1.0d-8).or.(ext_param%ktfin%value.gt.1.0d-8)) then
   N_mode=N_mode+1
endif

if (N_mode.eq.0) then
   stop 'ERROR: I have find no order parameter in the simulations'
else
   write(6,'(a,2x,I3,2x,a)') 'I have found',N_mode,'order parameters'
endif

get_num_mode=N_mode

end function get_num_mode

! routine that order the modes as fonction of what is already in the matrix
subroutine order_mode(i,name)
implicit none
integer, intent(in) :: i
character(len=*), intent(in) :: name
! internal

my_order_parameters(i)%name=name

if (i.eq.1) then
   write(6,'(2(a,2x))') 'first mode is',name
   my_order_parameters(i)%start=1
   my_order_parameters(i)%end=my_order_parameters(i)%start+dim_name_mode(name)-1

else
   my_order_parameters(i)%start=my_order_parameters(i-1)%end+1
   my_order_parameters(i)%end=my_order_parameters(i)%start+dim_name_mode(name)-1

endif

end subroutine order_mode



!!!!!!!!!!!!!!!!!!!!!!!
!
! function that gives if the mode has a dimension 1,2,3
!
!!!!!!!!!!!!!!!!!!!!!!!
integer function dim_name_mode(name)
implicit none
character(len=*), intent(in) :: name
dim_name_mode=0

select case (name)
  case ('magnetic')
    dim_name_mode=3
  case ('temperature')
    dim_name_mode=1
  case ('Bfield')
    dim_name_mode=3
  case ('Efield')
    dim_name_mode=3
  case ('polarisation')
    dim_name_mode=3
  case default
    stop 'this mode does not exist'
end select

end function dim_name_mode

end module m_lattice
