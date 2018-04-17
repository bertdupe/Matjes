module m_init_variables
use m_derived_types

private
public :: init_variables

interface init_variables
  module procedure init_simu
  module procedure init_lattice
  module procedure init_io_parameter
end interface

contains

!!!!!!!!!!!!!!!!!!
! initialization of the derived type variables
!!!!!!!!!!!!!!!!!!

! Initialization of the io parameters
subroutine init_io_parameter(io_simu)
implicit none
type(io_parameter), intent(out) :: io_simu

! io_of the simulation
io_simu%io_dispersion=.False.
io_simu%io_qorien=.False.
io_simu%io_spstmL=.False.
io_simu%io_spstmonly=.False.
io_simu%io_Xstruct=.False.
io_simu%io_fft_Xstruct=.False.
io_simu%io_topo=.False.
io_simu%io_topohall=.False.
io_simu%io_frequency=1
io_simu%io_warning=.True.
io_simu%io_writing=1

end subroutine init_io_parameter

! Initialization of the simulation type
subroutine init_simu(simu)
implicit none
type(bool_var), intent(out) :: simu
! internal variables

! type of the simulation
simu%value=.False.

! name of the simulation
simu%name=''

end subroutine init_simu

! Initialization of the lattice type

subroutine init_lattice(my_lattice)
implicit none
type(lattice), intent(out) :: my_lattice(:)
! internal
integer :: i,nsize
type(lattice) :: start_lattice

nsize=size(my_lattice)

start_lattice%areal=0.0
start_lattice%astar=0.0
start_lattice%alat=0.0

start_lattice%boundary=.True.

!nullify(start_lattice%l_modes)

do i=1,nsize
    my_lattice(i)=start_lattice
enddo

end subroutine init_lattice

end module m_init_variables
