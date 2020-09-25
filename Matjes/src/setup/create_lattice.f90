module m_lattice
use m_basic_types, only : var_name
use m_derived_types, only : order_parameter

! array that stores what are the order parameter used
type(order_parameter), protected, allocatable, save, target :: my_order_parameters(:)
! Array of spins on the lattice (matrix)

private
public :: create_lattice,my_order_parameters
contains

subroutine create_lattice(mag_lattice,motif,ext_param,nb_orbitals)
!subroutine that initializes the the lattice and its orderparameters
use m_derived_types
use m_type_lattice
implicit none
integer, intent(in) :: nb_orbitals
type(cell), intent(in) :: motif
type(lattice), intent(inout) :: mag_lattice
type(simulation_parameters), intent(in) :: ext_param
! internal
integer :: dim_lat(3), nmag, N_mode,dim_mode
! slopes
integer :: i,j,k,l

integer :: dim_mode_arr(number_different_order_parameters)

dim_lat=mag_lattice%dim_lat
nmag=1
mag_lattice%nmag=nmag
!find out how many order parameters are used
N_mode=get_num_mode(motif,ext_param,nb_orbitals)

!(old) initialize variable saving information what order parameter is what
allocate(my_order_parameters(N_mode))
do i=1,N_mode
    my_order_parameters(i)%name=''
    my_order_parameters(i)%start=-1
    my_order_parameters(i)%end=-1
enddo
!get sizes of the order parameters
Call get_dim_mode(motif,ext_param,nb_orbitals,dim_mode_arr)
dim_mode=sum(dim_mode_arr)

!old orderparameter format
mag_lattice%dim_mode=dim_mode
Call mag_lattice%ordpar%init(mag_lattice)

!new orderparameter format
if(dim_mode_arr(1)>0) Call mag_lattice%M%init(mag_lattice,dim_mode_arr(1))
if(dim_mode_arr(2)>0) Call mag_lattice%E%init(mag_lattice,dim_mode_arr(2))
if(dim_mode_arr(3)>0) Call mag_lattice%B%init(mag_lattice,dim_mode_arr(3))
if(dim_mode_arr(4)>0) Call mag_lattice%T%init(mag_lattice,dim_mode_arr(4))

end subroutine create_lattice

!!!!!!!!!!!!!!!!!!!!!!!!
! function that gets the number of order parameter
!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_dim_mode(motif,ext_param,nb_orbitals,dim_mode_arr)
use m_derived_types
implicit none
integer, intent(in) :: nb_orbitals
type(cell), intent(in) :: motif
type(simulation_parameters), intent(in) :: ext_param
integer :: dim_mode_arr(number_different_order_parameters)
! internal parameters
integer :: nmag,N_mode
real(kind=8) :: Field(3)

dim_mode_arr=0
nmag=count(motif%atomic(:)%moment.gt.0.0d0)
dim_mode_arr(1)=nmag*3

! for each order parameter, organize the matrix of the mode
N_mode=0
if (nmag.ne.0) then
    N_mode=N_mode+1
    call order_mode(N_mode,'magnetic',3*nmag)
endif

! electric field
Field=ext_param%E_ext%value
if (sqrt(Field(1)**2+Field(2)**2+Field(3)**2).gt.1.0d-8) then
   dim_mode_arr(2)=3
   N_mode=N_mode+1
   call order_mode(N_mode,'Efield',3)
   write(6,'(a)') 'electric field found'
   write(6,'(a,2x,I5)') 'Hamiltonian has the dimension',sum(dim_mode_arr)**2
endif

! Magnetic field
Field=ext_param%H_ext%value
if (sqrt(Field(1)**2+Field(2)**2+Field(3)**2).gt.1.0d-8) then
   dim_mode_arr(3)=3
   N_mode=N_mode+1
   call order_mode(N_mode,'Bfield',3)
   write(6,'(a)') 'magnetic field found'
   write(6,'(a,2x,I5)') 'Hamiltonian has the dimension',sum(dim_mode_arr)**2
endif

! Temperature
if ((ext_param%ktini%value.gt.1.0d-8).or.(ext_param%ktfin%value.gt.1.0d-8)) then
   dim_mode_arr(4)=1
   N_mode=N_mode+1
   call order_mode(N_mode,'temperature',1)
   write(6,'(a)') 'Temperature found'
   write(6,'(a,2x,I5)') 'Hamiltonian has the dimension',sum(dim_mode_arr)**2
endif

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!
! function that gets the number of different order parameter
!!!!!!!!!!!!!!!!!!!!!!!!
integer function get_num_mode(motif,ext_param,nb_orbitals)
use m_derived_types
implicit none
integer, intent(in) :: nb_orbitals
type(cell), intent(in) :: motif
type(simulation_parameters), intent(in) :: ext_param
! internal parameters
integer :: N_mode,nmag
real(kind=8) :: Field(3)

N_mode=0
! counts one mode per magnetic atom in the unit cell
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

! routine that order the modes as function of what is already in the matrix
subroutine order_mode(i,name,dim_mode)
implicit none
integer, intent(in) :: i, dim_mode
character(len=*), intent(in) :: name
! internal

my_order_parameters(i)%name=name

if (i.eq.1) then
   write(6,'(2(a,2x))') 'first mode is',name
   my_order_parameters(i)%start=1
   my_order_parameters(i)%end=my_order_parameters(i)%start+dim_mode-1
else
   my_order_parameters(i)%start=my_order_parameters(i-1)%end+1
   my_order_parameters(i)%end=my_order_parameters(i)%start+dim_mode-1

endif

end subroutine order_mode

end module m_lattice
