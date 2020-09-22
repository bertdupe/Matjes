module m_init_heavyside
use m_derived_types
implicit none
private
public :: init_heavyside
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the starting configuration as a heavyside function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_heavyside(my_lattice,my_motif,start,end)
type (lattice), intent(inout) :: my_lattice
type (cell), intent(in) :: my_motif
integer, intent(in) :: start,end
! internal variable
integer :: i_z,i_y,i_x,i_m,Nx,Ny,Nz,size_mag

! get the position of the sites on the lattice
Nx=my_lattice%dim_lat(1)
Ny=my_lattice%dim_lat(2)
Nz=my_lattice%dim_lat(3)
size_mag=count(my_motif%atomic(:)%moment.gt.0.0d0)


do i_m=1,size_mag
   do i_z=1,Nz
      do i_y=1,Ny
         do i_x=1,Nx

            my_lattice%ordpar%l_modes(i_x,i_y,i_z,i_m)%w(start:start+1)=0.0d0
            my_lattice%ordpar%l_modes(i_x,i_y,i_z,i_m)%w(end)=1-2*(2*i_x/Nx)

          enddo
      enddo
   enddo
enddo

end subroutine init_heavyside

end module m_init_heavyside
