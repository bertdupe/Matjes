module m_init_heavyside
use m_derived_types
implicit none
private
public :: init_heavyside
interface init_heavyside
    module procedure init_heavyside_old
    module procedure init_heavyside_new
end interface
contains


subroutine init_heavyside_new(lat,dim_mode,state)
    type (lattice),intent(in)       :: lat
    integer, intent(in)             :: dim_mode
    real(8),pointer,intent(inout)   :: state(:)
    ! internal variable
    integer         :: nmag
    real(8),pointer :: state_x(:,:,:,:)
   
    nmag=lat%cell%num_mag()
    state_x(1:3,1:nmag,1:lat%dim_lat(1),1:lat%dim_lat(2)*lat%dim_lat(3))=>state

    state_x=0.0d0
    state_x(3,:,:,:)=-1.0d0
    state_x(3,:,1:lat%dim_lat(1)/2,:)=1.0d0
end subroutine 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the starting configuration as a heavyside function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_heavyside_old(my_lattice,my_motif,start,end)
type (lattice), intent(inout) :: my_lattice
type(t_cell), intent(in) :: my_motif
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

end subroutine 

end module m_init_heavyside
