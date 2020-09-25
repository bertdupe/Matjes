module m_init_DW
use m_derived_types
implicit none
private
public :: init_DW

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the starting configuration as a domain wall along the x direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_DW(my_lattice,my_motif,start,end)
type (lattice), intent(inout) :: my_lattice
type(t_cell), intent(in) :: my_motif
integer, intent(in) ::  start,end
! internal variables
integer :: i_z,i_y,i_x,i_m,Nx,Ny,Nz,size_mag,dw_position,i_w,shape_lattice(4)
real(kind=8) :: alpha

! get the position of the sites on the lattice
Nx=my_lattice%dim_lat(1)
Ny=my_lattice%dim_lat(2)
Nz=my_lattice%dim_lat(3)
shape_lattice=shape(my_lattice%ordpar%l_modes)
size_mag=shape_lattice(4)

do i_x=1,Nx
   if ( (2*i_x/Nx) == 1) then
      dw_position=i_x
      exit
   endif
enddo

do i_x=1,Nx
   do i_w=-4,4
      alpha=real(5+i_w)/10.0d0*acos(-1.0d0)
      if ( i_x+i_w == dw_position ) then
         do i_y=1,Ny
            do i_z=1,Nz
               do i_m=1,size_mag

           my_lattice%ordpar%l_modes(i_x,i_y,i_z,i_m)%w(start)=sin(alpha)
           my_lattice%ordpar%l_modes(i_x,i_y,i_z,i_m)%w(start+1)=0.0d0
           my_lattice%ordpar%l_modes(i_x,i_y,i_z,i_m)%w(end)=-1.0d0*cos(alpha)

               enddo
            enddo
         enddo

      elseif ( i_x+i_w .gt. dw_position ) then
         do i_y=1,Ny
            do i_z=1,Nz
               do i_m=1,size_mag

           my_lattice%ordpar%l_modes(i_x,i_y,i_z,i_m)%w(start:start+1)=0.0d0
           my_lattice%ordpar%l_modes(i_x,i_y,i_z,i_m)%w(end)=-1.0d0

               enddo
            enddo
         enddo

      else
         do i_y=1,Ny
            do i_z=1,Nz
               do i_m=1,size_mag

           my_lattice%ordpar%l_modes(i_x,i_y,i_z,i_m)%w(start:start+1)=0.0d0
           my_lattice%ordpar%l_modes(i_x,i_y,i_z,i_m)%w(end)=1.0d0

               enddo
            enddo
         enddo
      endif
   enddo
enddo

end subroutine init_DW

end module m_init_DW
