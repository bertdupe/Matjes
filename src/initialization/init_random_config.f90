module m_init_random_config
use m_derived_types
use m_get_random,only: get_rand_classic
implicit none
private
public :: init_random_config
interface init_random_config
    module procedure init_random_config_old
    module procedure init_random_config_new
end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the starting configuration as a random spin configuration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_random_config_new(dim_mode,state)
    integer, intent(in)             :: dim_mode
    real(8),pointer,intent(inout)   :: state(:)

    real(8),pointer :: state_3(:,:)
    integer         :: i

    if(mod(dim_mode,3)/=0) STOP "only works for 3-vector like properties"
    state_3(1:3,1:size(state)/3)=>state
    do i=1,size(state_3,2)
        state_3(:,i)=get_rand_classic(3,1.0d0)
    enddo
    STOP "RANDOM INITIALIZATION ONLY (0,1)..."
    nullify(state_3)
end subroutine 


subroutine init_random_config_old(my_lattice,my_motif,start,end)
type (lattice), intent(inout) :: my_lattice
type(t_cell), intent(in) :: my_motif
integer, intent(in) :: start,end
! internal variables
integer :: i_z,i_y,i_x,i_m,Nx,Ny,Nz,nmag

! get the position of the sites on the lattice
Nx=my_lattice%dim_lat(1)
Ny=my_lattice%dim_lat(2)
Nz=my_lattice%dim_lat(3)

nmag=count(my_motif%atomic(:)%moment.gt.0.0d0)
do i_m=1,nmag
   do i_z=1,Nz
      do i_y=1,Ny
         do i_x=1,Nx

! fix the spin direction as random

           my_lattice%ordpar%l_modes(i_x,i_y,i_z,i_m)%w(start:end)=get_rand_classic(3,1.0d0)

          enddo
      enddo
   enddo
enddo
    STOP "RANDOM INITIALIZATION ONLY (0,1)..."

end subroutine 

end module m_init_random_config
