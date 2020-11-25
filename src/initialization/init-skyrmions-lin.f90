module m_init_sky_lin
use m_derived_types
implicit none
private
public :: init_sky_lin
interface init_sky_lin
    module procedure init_sky_lin_old
    module procedure init_sky_lin_new
end interface
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sets the magnetization according the skyrmion structure of:
!   Garnier, M., Mesaros, A., & Simon, P. (2019). 
!   Topological superconductivity with deformable magnetic skyrmions.
!   Communications Physics, 2(1), 126.
!   https://doi.org/10.1038/s42005-019-0226-5
! i.e.
! relative to some center we sets
! n(r)=[sin(f(r))*cos(q*phi),sin(f(r))*sin(q*phi),cos(f(r))]
! with f(r) being a linear function with a given slope(multplied with 2*pi from the input),
! q being an integer giving the azimuthal winding number,
! and phi being the polar angle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_sky_lin_new(io,fname,lat,ordname,dim_mode,state)
    use m_io_utils, only: get_parameter
    use m_init_util, only: get_pos_vec
    use m_constants, only : pi
    integer,intent(in)              :: io       !init-file io-unit
    character(*),intent(in)         :: fname    !init-file name 
    type(lattice), intent(in)       :: lat      !entire lattice containing geometric information
    character(*),intent(in)         :: ordname  !name of the order parameter
    integer,intent(in)              :: dim_mode !dimension of the order parameter in each cell
    real(8),pointer,intent(inout)   :: state(:) !pointer the the order parameter

    real(8),allocatable,target :: pos(:)
    real(8),pointer :: pos_3(:,:),state_3(:,:)
    real(8)     ::  f_r,qphi
    integer     :: i

    real(8)     ::  slope,center(3)
    integer     ::  q

    !get parameters for punch
    center=0.0d0
    q=1
    slope=10.0
    !older input without ordername? (legacy?)
    call get_parameter(io,fname,'skylin_pos',3,center)
    call get_parameter(io,fname,'skylin_slope',slope)
    call get_parameter(io,fname,'skylin_q',q)
    call get_parameter(io,fname,'skylin_pos_'//ordname,3,center)
    call get_parameter(io,fname,'skylin_slope_'//ordname,slope)
    call get_parameter(io,fname,'skylin_q_'//ordname,q)
    slope=2.0d0*pi/slope
    
    Call get_pos_vec(lat,dim_mode,ordname,pos)
    pos_3(1:3,1:size(pos)/3)=>pos
    state_3(1:3,1:size(pos)/3)=>state

    do i=1,size(pos_3,2)
        pos_3(:,i)=pos_3(:,i)-center
    enddo
    
    do i=1,size(pos_3,2)
        f_r=norm2(pos_3(:,i))*slope
        qphi=atan2(pos_3(2,i),pos_3(1,i))*real(q,8)
        state_3(1,i)= sin(f_r)*cos(qphi)
        state_3(2,i)= sin(f_r)*sin(qphi)
        state_3(3,i)= cos(f_r)
    enddo

    nullify(pos_3,state_3)
    deallocate(pos)
end subroutine 


subroutine init_sky_lin_old(io,fname,my_lattice,my_motif,m_start,m_end)
use m_io_utils
type (lattice), intent(inout) :: my_lattice
type(t_cell), intent(in) :: my_motif
integer, intent(in) :: io,m_start,m_end
character(len=*), intent(in) :: fname
! internal variables
real(8)     ::  pos(3)
real(8)     ::  slope
integer     ::  q
real(8),parameter   ::  pi=3.1415926535897932384626433d0
character(len=30) :: configuration

!get parameters for punch
pos=0.0d0
q=1
slope=10.0
call get_parameter(io,fname,'skylin_pos',3,pos)
call get_parameter(io,fname,'skylin_slope',slope)
call get_parameter(io,fname,'skylin_q',q)
slope=2.0d0*pi/slope

Call set_sky_lin(pos,slope,q,my_lattice,my_motif,m_start,m_end)

end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!sets all magnetization to zero that is farther from "center" than radius
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_sky_lin(center,slope,q,my_lattice,my_motif,m_start,m_end)
use m_get_position, only: get_position
Implicit None
type (lattice), intent(inout) :: my_lattice
type(t_cell), intent(in) :: my_motif
! center of circle region
real(kind=8), intent(in) :: center(3)
! slope for the radial part
real(kind=8),intent(in) :: slope
!azimuthal winding number
integer,intent(in)      :: q

integer, intent(in) :: m_start,m_end
! Internal variables
Integer:: i_x,i_y,i_z,i_m,dim_lat(3),nmag
Integer:: i
! get the position of the sites
real(kind=8), allocatable,target :: pos(:,:,:,:,:)
real(8),pointer         ::  pos_mat2(:,:)
integer :: Nx,Ny,Nz

real(8)     ::  f_r,qphi

dim_lat=my_lattice%dim_lat
nmag=count(my_motif%atomic(:)%moment.gt.0.0d0)
Nx=dim_lat(1); Ny=dim_lat(2); Nz=dim_lat(3)

!get all positions("pos") relative to center
allocate(pos(3,Nx,Ny,Nz,nmag),source=0.0d0)
call get_position(pos,dim_lat,my_lattice%areal,my_motif)
pos_mat2(1:3,1:Nx*Ny*Nz*nmag)=>pos
do i=1,size(pos_mat2,2)
    pos_mat2(:,i)=pos_mat2(:,i)-center
enddo
nullify(pos_mat2)

do i_m=1,nmag
   Do i_z=1,dim_lat(3)
      Do i_y=1,dim_lat(2)
         Do i_x=1,dim_lat(1)
           f_r=norm2(pos(:,i_x,i_y,i_z,i_m))*slope
           qphi=atan2(pos(2,i_x,i_y,i_z,i_m),pos(1,i_x,i_y,i_z,i_m))*real(q,8)
           my_lattice%ordpar%l_modes(i_x,i_y,i_z,i_m)%w(m_start)   = sin(f_r)*cos(qphi)
           my_lattice%ordpar%l_modes(i_x,i_y,i_z,i_m)%w(m_start+1) = sin(f_r)*sin(qphi)
           my_lattice%ordpar%l_modes(i_x,i_y,i_z,i_m)%w(m_start+2) = cos(f_r)
         enddo
      enddo
   enddo
enddo
end subroutine

end module 
