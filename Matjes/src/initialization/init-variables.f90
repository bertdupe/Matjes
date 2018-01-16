module m_init_variables
use m_derived_types

interface init_variables
  module procedure init_simu
  module procedure init_lattice
end interface

private
public :: init_variables,init_lattice

contains

!!!!!!!!!!!!!!!!!!
! initialization of the derived type variables
!!!!!!!!!!!!!!!!!!

! Initialization of the simulation type
subroutine init_simu(simu)
implicit none
type(bool_var), intent(out) :: simu(:)
! internal variables
integer :: i, n_size

n_size=size(simu)

do i=1,n_size
! type of the simulation
   simu(i)%value=.False.

! name of the simulations
   simu(i)%name=''

! name of the variables
   simu(i)%var_name=''

enddo

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

nullify(start_lattice%l_modes)

do i=1,nsize
    my_lattice(i)=start_lattice
enddo

end subroutine init_lattice

end module m_init_variables
