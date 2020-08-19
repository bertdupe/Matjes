module m_init_punch
use m_derived_types
implicit none
private
public :: init_punch
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Sets magnetization to zero outside if inserted region
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_punch(io,fname,my_lattice,my_motif,m_start,m_end)
!punches out an area with magnetization (i.e. removes magnetization outsize of the region)
use m_io_utils
type (lattice), intent(inout) :: my_lattice
type (cell), intent(in) :: my_motif
integer, intent(in) :: io,m_start,m_end
character(len=*), intent(in) :: fname
! internal variables
real(8)     ::  pos(3),radius

!get parameters for punch
pos=0.0d0
radius=Huge(1.0d0)

call get_parameter(io,fname,'punch_posx',pos(1))
call get_parameter(io,fname,'punch_posy',pos(2))
call get_parameter(io,fname,'punch_posz',pos(3))
call get_parameter(io,fname,'punch_radius',radius)

call punch_texture_circle(pos,radius,my_lattice,my_motif,m_start,m_end)

end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!sets all magnetization to zero that is farther from "center" than radius
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine punch_texture_circle(center,radius,my_lattice,my_motif,m_start,m_end)
use m_get_position, only: get_position
Implicit None
type (lattice), intent(inout) :: my_lattice
type (cell), intent(in) :: my_motif
! center of circle region
real(kind=8), intent(in) :: center(3)
! radius of circular region
real(kind=8),intent(in) :: radius
integer, intent(in) :: m_start,m_end
! Internal variables
Integer:: i_x,i_y,i_z,i_m,dim_lat(3),nmag
Integer:: i
! get the position of the sites
real(kind=8), allocatable,target :: pos(:,:,:,:,:)
real(8),pointer         ::  pos_mat2(:,:)
integer :: Nx,Ny,Nz

dim_lat=my_lattice%dim_lat
nmag=count(my_motif%atomic(:)%moment.gt.0.0d0)
Nx=dim_lat(1); Ny=dim_lat(2); Nz=dim_lat(3)

!get all positions("pos") relative to center
allocate(pos(3,Nx,Ny,Nz,nmag),source=0.0d0)
call get_position(pos,dim_lat,my_lattice%areal,my_motif)
pos_mat2(1:3,1:Nx*Ny*Nz*nmag)=>pos
pos_mat2=matmul(my_lattice%areal,pos_mat2)
do i=1,size(pos_mat2,2)
    pos_mat2(:,i)=pos_mat2(:,i)-center
enddo
nullify(pos_mat2)

do i_m=1,nmag
   Do i_z=1,dim_lat(3)
      Do i_y=1,dim_lat(2)
         Do i_x=1,dim_lat(1)
            if(norm2(pos(:,i_x,i_y,i_z,i_m))>radius)then
               my_lattice%l_modes(i_x,i_y,i_z,i_m)%w(m_start:m_end) = 0.0d0
            endif
         enddo
      enddo
   enddo
enddo
end subroutine

end module 
