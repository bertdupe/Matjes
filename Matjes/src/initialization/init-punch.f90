module m_init_punch
use m_derived_types
implicit none
private
public :: init_punch
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Sets magnetization to zero outside if inserted region
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_punch(io,fname,my_lattice,my_motif,mode_name,m_start,m_end)
!punches out an area with magnetization (i.e. removes magnetization outsize of the region)
use m_io_utils
type (lattice), intent(inout) :: my_lattice
type (cell), intent(in) :: my_motif
integer, intent(in) :: io,m_start,m_end
character(len=*), intent(in) :: fname,mode_name
! internal variables
real(8)     ::  pos(3),radius,width,height,angle
character(len=30) :: configuration

select case(adjustl(mode_name))
  case('magnetic')
    continue
  case default
    return
end select

!get parameters for punch
pos=0.0d0
radius=Huge(1.0d0)
width=Huge(1.0d0)
height=Huge(1.0d0)
angle=0.0d0
configuration="nothing"
call get_parameter(io,fname,'punch_config',configuration)
call get_parameter(io,fname,'punch_pos',3,pos)
call get_parameter(io,fname,'punch_radius',radius)
call get_parameter(io,fname,'punch_width',width)
call get_parameter(io,fname,'punch_height',height)
call get_parameter(io,fname,'punch_angle',angle)

select case (adjustl(configuration))
  case('circle')
    call punch_texture_circle(pos,radius,my_lattice,my_motif,m_start,m_end)
  case('hexagon')
    Call punch_texture_hexagon(pos,width,angle,my_lattice,my_motif,m_start,m_end)
  case('rectangle')
    Call punch_texture_rectangle(pos,width,height,angle,my_lattice,my_motif,m_start,m_end)
end select


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
integer                 ::  Nx,Ny,Nz
integer                 ::  Npos

dim_lat=my_lattice%dim_lat
nmag=count(my_motif%atomic(:)%moment.gt.0.0d0)
Nx=dim_lat(1); Ny=dim_lat(2); Nz=dim_lat(3)

!get all positions("pos") relative to center
allocate(pos(3,Nx,Ny,Nz,nmag),source=0.0d0)
call get_position(pos,dim_lat,my_lattice%areal,my_motif)
Npos=Nx*Ny*Nz*nmag
pos_mat2(1:3,1:Npos)=>pos
do i=1,Npos
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


subroutine punch_texture_rectangle(center,width,height,angle_in,my_lattice,my_motif,m_start,m_end)
use m_get_position, only: get_position
Implicit None
type (lattice), intent(inout) :: my_lattice
type (cell), intent(in) :: my_motif
!rectangle parameters
real(kind=8), intent(in) :: center(3)
real(kind=8),intent(in) :: width,height,angle_in
integer, intent(in) :: m_start,m_end

! Internal variables
Integer:: i_x,i_y,i_z,i_m,dim_lat(3),nmag
Integer:: i
! get the position of the sites
real(kind=8), allocatable,target :: pos(:,:,:,:,:)
real(8),pointer         ::  pos_mat2(:,:)
integer                 ::  Nx,Ny,Nz
integer                 ::  Npos
!rotation parameters
real(8)                 ::  angle
real(8)                 ::  rot(3,3)
real(8),parameter       ::  pi=3.1415926535897932384626433d0
dim_lat=my_lattice%dim_lat
nmag=count(my_motif%atomic(:)%moment.gt.0.0d0)
Nx=dim_lat(1); Ny=dim_lat(2); Nz=dim_lat(3)

!get all positions("pos") relative to center
allocate(pos(3,Nx,Ny,Nz,nmag),source=0.0d0)
call get_position(pos,dim_lat,my_lattice%areal,my_motif)
Npos=Nx*Ny*Nz*nmag
pos_mat2(1:3,1:Npos)=>pos
do i=1,Npos
    pos_mat2(:,i)=pos_mat2(:,i)-center
enddo

!rotate position
angle=angle_in*2.0d0*pi/360.0
rot(:,1)=[ cos(angle),-sin(angle),0.0d0]
rot(:,2)=[ sin(angle), cos(angle),0.0d0]
rot(:,3)=[0.0d0,0.0d0,1.0d0]
pos_mat2=matmul(rot,pos_mat2)

nullify(pos_mat2)

do i_m=1,nmag
   Do i_z=1,dim_lat(3)
      Do i_y=1,dim_lat(2)
         Do i_x=1,dim_lat(1)
            if(abs(pos(1,i_x,i_y,i_z,i_m))>width*0.5d0.or.abs(pos(2,i_x,i_y,i_z,i_m))>height*0.5d0)then
               my_lattice%l_modes(i_x,i_y,i_z,i_m)%w(m_start:m_end) = 0.0d0
            endif
         enddo
      enddo
   enddo
enddo
end subroutine


subroutine punch_texture_hexagon(center,width,angle_in,my_lattice,my_motif,m_start,m_end)
use m_get_position, only: get_position
Implicit None
type (lattice), intent(inout) :: my_lattice
type (cell), intent(in) :: my_motif
! center of circle region
real(kind=8), intent(in) :: center(3)
! radius of circular region
real(kind=8),intent(in) :: width,angle_in !minor width 
integer, intent(in) :: m_start,m_end
! Internal variables
Integer:: i_x,i_y,i_z,i_m,dim_lat(3),nmag
Integer:: i
! get the position of the sites
real(kind=8), allocatable,target :: pos(:,:,:,:,:)
real(8),pointer         ::  pos_mat2(:,:)
integer :: Nx,Ny,Nz

integer                 ::  Npos
logical,allocatable,target     ::  outside(:)
logical,pointer         ::  outside_point(:,:,:,:)

real(8)                 ::  angle
real(8)                 ::  rot(3,3)
real(8),parameter       ::  pi=3.1415926535897932384626433d0

dim_lat=my_lattice%dim_lat
nmag=count(my_motif%atomic(:)%moment.gt.0.0d0)
Nx=dim_lat(1); Ny=dim_lat(2); Nz=dim_lat(3)

!get all positions("pos") relative to center
allocate(pos(3,Nx,Ny,Nz,nmag),source=0.0d0)
Npos=Nx*Ny*Nz*nmag
call get_position(pos,dim_lat,my_lattice%areal,my_motif)
pos_mat2(1:3,1:Npos)=>pos
do i=1,Npos
    pos_mat2(:,i)=pos_mat2(:,i)-center
enddo

!rotate position (e.g. 30 deg for other orientation) 
angle=angle_in*2.0d0*pi/360.0
rot(:,1)=[ cos(angle),-sin(angle),0.0d0]
rot(:,2)=[ sin(angle), cos(angle),0.0d0]
rot(:,3)=[0.0d0,0.0d0,1.0d0]
pos_mat2=matmul(rot,pos_mat2)

!only check upper quadrant
pos_mat2=abs(pos_mat2)/width 
allocate(outside(Npos),source=.False.)
!check if outside maximal rectangle
outside=pos_mat2(1,:)>=sqrt(3.0d0)*0.5d0.or.pos_mat2(2,:)>=1.0d0
!check is below 
outside=outside.or.pos_mat2(2,:)-1.0d0+pos_mat2(1,:)/sqrt(3.0d0)>0.0d0
nullify(pos_mat2)

!set all outside==.true. magnetizations to zero
outside_point(1:Nx,1:Ny,1:Nz,1:nmag)=>outside
do i_m=1,nmag
   Do i_z=1,dim_lat(3)
      Do i_y=1,dim_lat(2)
         Do i_x=1,dim_lat(1)
            if(outside_point(i_x,i_y,i_z,i_m))then
               my_lattice%l_modes(i_x,i_y,i_z,i_m)%w(m_start:m_end) = 0.0d0
            endif
         enddo
      enddo
   enddo
enddo
nullify(outside_point)
end subroutine


end module 
