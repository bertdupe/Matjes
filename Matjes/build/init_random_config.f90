module m_init_random_config
use m_derived_types
use m_get_random
use mtprng, only : mtprng_state
implicit none
private
public :: init_random_config

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the starting configuration as a random spin configuration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_random_config(my_lattice,my_motif,state)
type (lattice), intent(inout) :: my_lattice
type (cell), intent(in) :: my_motif
type (mtprng_state),intent(inout) :: state
! internal variables
integer :: i_z,i_y,i_x,i_m,Nx,Ny,Nz,size_mag,nmag

! get the position of the sites on the lattice
Nx=my_lattice%dim_lat(1)
Ny=my_lattice%dim_lat(2)
Nz=my_lattice%dim_lat(3)
size_mag=size(my_motif%i_mom)

nmag=count(my_motif%atomic(:)%moment.gt.0.0d0)
do i_m=1,nmag
   do i_z=1,Nz
      do i_y=1,Ny
         do i_x=1,Nx

! fix the spin direction as random

           my_lattice%l_modes(i_x,i_y,i_z,i_m)%w(:)=get_rand_classic(state,3,1.0d0)

          enddo
      enddo
   enddo
enddo

end subroutine init_random_config

end module m_init_random_config
