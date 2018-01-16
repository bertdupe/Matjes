module m_init_modes
use m_derived_types
interface get_init_modes
module procedure init_3Dmodes
end interface get_init_modes

private
public :: get_init_modes

contains

subroutine init_3Dmodes(fname,my_lattice,motif)
use m_io_utils
implicit none
type (lattice), intent(inout) :: my_lattice
type (cell), intent(in) :: motif
character(len=*), intent(in) :: fname
! internal variables
integer :: dim_lat(3),n_column,io,nmag
! slope variables
integer :: i_x,i_y,i_z,i_m,j_lat

dim_lat=my_lattice%dim_lat
n_column=get_cols(fname)
nmag=count(motif%i_mom)

open(newunit=io,file=fname,form='formatted',status='old',action='read')
rewind(io)
!cc check for old format

select case(n_column)

case(3)

do i_m=1,nmag
   do i_z=1,dim_lat(3)
      do i_y=1,dim_lat(2)
         do i_x=1,dim_lat(1)
           read(io,*) (my_lattice%l_modes(i_x,i_y,i_z,i_m)%w(j_lat),j_lat=1,3)
         enddo
      enddo
   enddo
enddo

case(5)

do i_m=1,nmag
   do i_z=1,dim_lat(3)
      do i_y=1,dim_lat(2)
         do i_x=1,dim_lat(1)
           read(io,*) (my_lattice%l_modes(i_x,i_y,i_z,i_m)%w(j_lat),j_lat=3,5)
         enddo
      enddo
   enddo
enddo

case(6)

do i_m=1,nmag
   do i_z=1,dim_lat(3)
      do i_y=1,dim_lat(2)
         do i_x=1,dim_lat(1)
           read(io,*) (my_lattice%l_modes(i_x,i_y,i_z,i_m)%w(j_lat),j_lat=4,6)
         enddo
      enddo
   enddo
enddo

end select

close(io)

end subroutine init_3Dmodes

end module m_init_modes
