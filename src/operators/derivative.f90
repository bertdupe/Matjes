module m_derivative
use m_basic_types, only : vec_point
use m_derived_types, only : lattice


type(vec_point), allocatable :: derivative(:,:,:)

interface get_derivative
  module procedure get_derivative_int,get_derivative_out
end interface

interface calculate_derivative
  module procedure calculate_derivative_vector_field,calculate_derivative_scalar_field
end interface

private
public :: get_derivative,calculate_derivative
contains


!
! prepare the derivative in the module
!
subroutine get_derivative_int(field,mag_lattice)
use m_io_files_utils
use m_io_utils
implicit none
type(vec_point), target, intent(in) :: field(:)
type(lattice), intent(in) :: mag_lattice
! internal
integer, allocatable :: pos(:,:,:)
integer :: N_cell,i,n,dim_lat(4),N_world,j,k,location
integer :: io_input
logical :: periodic(3)
real(kind=8) :: DF(3,2)

dim_lat(1:3)=mag_lattice%dim_lat
dim_lat(4)=1
periodic=mag_lattice%periodic
N_world=size(mag_lattice%world)

n=1
io_input=open_file_read('input')
call get_parameter(io_input,'input','order_derivative',n)
call close_file('input',io_input)

N_cell=size(field)

allocate(pos(2*n+1,3,N_cell))
pos=0
call get_voisin(pos,dim_lat,periodic)

allocate(derivative(2*n+1,N_world,N_cell))

do i=1,N_cell
   do j=1,N_world
      do k=1,2*n+1
        location=pos(k,j,i)
        derivative(k,j,i)%w=>field(location)%w
      enddo
   enddo
enddo

io_input=open_file_write('derivative.dat')
do i=1,N_cell
  DF=0.0d0
  call calculate_derivative(DF,i)
enddo
call close_file('derivative.dat',io_input)

end subroutine get_derivative_int


!
! prepare the derivative and send it out of the module
!

subroutine get_derivative_out(field,mag_lattice,derivat_point)
implicit none
real(kind=8), target, intent(in) :: field(:,:)
type(lattice), intent(in) :: mag_lattice
type(vec_point), intent(inout) :: derivat_point(:,:,:)
! internal
integer, allocatable :: pos(:,:,:)
integer :: N_cell,i,dim_lat(4),N_world,j,k,location,shape_D(3),N_vois
logical :: periodic(3)

dim_lat(1:3)=mag_lattice%dim_lat
dim_lat(4)=1
periodic=mag_lattice%periodic
N_world=size(mag_lattice%world)
shape_D=shape(derivat_point)
N_cell=shape_D(3)
N_vois=shape_D(1)

allocate(pos(N_vois,3,N_cell))
pos=0
call get_voisin(pos,dim_lat,periodic)

do i=1,N_cell
   do j=1,N_world
      do k=1,N_vois
        location=pos(k,j,i)
        derivat_point(k,j,i)%w=>field(:,location)
      enddo
   enddo
enddo

end subroutine get_derivative_out

!
! derivative of a vector field
!
!
subroutine calculate_derivative_vector_field(DF,iomp)
use m_vector, only : cross,norm
implicit none
integer, intent(in) :: iomp
real(kind=8), intent(inout) :: DF(:,:)
! internal
integer :: shape_field(2),i,j
real(kind=8) :: step(3),norm_int,step_DF(3)

shape_field=shape(derivative(:,:,iomp))

do i=1,shape_field(2)
   do j=1,shape_field(1)-1

      step=cross(derivative(j,i,iomp)%w,derivative(j+1,i,iomp)%w)
      norm_int=norm(step)
      if (norm_int.lt.1.0d-7) then
         DF(:,i)=0.0
         cycle
      endif
      step=step/norm_int
      step_DF=-cross(derivative(j,i,iomp)%w,step)
      norm_int=norm(step_DF)
      DF(:,i)=step_DF/norm_int

   enddo

   DF(:,i)=DF(:,i)/real(shape_field(1)-1)
enddo

end subroutine calculate_derivative_vector_field

!
! derivative of a scalar field
!
!
subroutine calculate_derivative_scalar_field(DF,field,n_dim,n_vois)
implicit none
integer, intent(in) :: n_dim,n_vois
type(vec_point), intent(in) :: field(:,:)
real(kind=8), intent(inout) :: DF(:,:)
! internal
integer :: i,j

DF=0.0d0

do i=1,n_dim
  do j=1,n_vois-1

     DF(:,i)=DF(:,i)+field(j+1,i)%w-field(j,i)%w

   enddo
enddo

DF=DF/real(n_vois-1)

end subroutine calculate_derivative_scalar_field





!
! get the neihbors in the direction u, v and w
!
!
subroutine get_voisin(pos,dim_lat,periodic)
use m_get_position
implicit none
integer, intent(inout) :: pos(:,:,:)
integer, intent(in) :: dim_lat(:)
logical, intent(in) :: periodic(:)
! internal
integer :: i_x,i_y,i_z,i_m,N_vois,j
integer :: ipu,ipv,ipw,Ilat(4),i_pos,i_pos_v

N_vois=(size(pos,1)-1)/2

do i_m=1,dim_lat(4)
   do i_z=1,dim_lat(3)
      do i_y=1,dim_lat(2)
         do i_x=1,dim_lat(1)
            Ilat=(/i_x,i_y,i_z,i_m/)
            i_pos=get_position_ND_to_1D(Ilat,dim_lat)

            do j=-N_vois,N_vois

               ipu=neighbor(i_x,j,Periodic(1),dim_lat(1))
               Ilat(1)=ipu
               i_pos_v=get_position_ND_to_1D(Ilat,dim_lat)
               pos(j+N_vois+1,1,i_pos)=i_pos_v

               ipv=neighbor(i_y,j,Periodic(2),dim_lat(2))
               Ilat(1)=i_x
               Ilat(2)=ipv
               i_pos_v=get_position_ND_to_1D(Ilat,dim_lat)
               pos(j+N_vois+1,2,i_pos)=i_pos_v

               ipw=neighbor(i_z,j,Periodic(3),dim_lat(3))
               Ilat(2)=i_y
               Ilat(3)=ipw
               i_pos_v=get_position_ND_to_1D(Ilat,dim_lat)
               pos(j+N_vois+1,3,i_pos)=i_pos_v
            enddo

         enddo
      enddo
   enddo
enddo

end subroutine

integer function neighbor(i,t,Periodic,N_max)
implicit none
integer :: i,N_max,t
logical :: Periodic

if (Periodic) then
   neighbor=mod(i+t-1+N_max,N_max)+1
elseif ((i+t.gt.N_max).or.(i+t.lt.1)) then
   neighbor=i
else
   neighbor=i+t
endif

end function

end module m_derivative
