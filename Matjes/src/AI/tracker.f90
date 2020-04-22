module m_tracker
use m_topoplot
use m_topo_commons
use m_derived_types, only : lattice
use m_basic_types, only : vec_point
use m_derivative

! tracking parameters
type pattern
  integer, allocatable :: interior(:)
  integer, allocatable :: boundary(:)
  integer :: center,Nsites
end type

! derivation of the field of interest
type(vec_point), allocatable :: derivative(:,:,:)

! the field of interest
real(kind=8), allocatable :: field(:,:)

type(pattern), allocatable :: all_patterns(:)

! tracking frequency
integer :: track_freq=1

private
public :: init_tracking, plot_tracking

contains

!
! First guess of what should be tracked
!
subroutine init_tracking(my_lattice)
use m_constants, only : pi
use m_io_files_utils
use m_io_utils
implicit none
type(lattice), intent(in) :: my_lattice
! internal variable
integer :: i,N_pat,Nsites,N_cell,N_world,dim_lat(3),io_input
! derivative parameter
integer :: order_derivative
real(kind=8) :: Q(5)

write(6,'(a)') ''
write(6,'(a)') 'initialization of the tracking procedure'
write(6,'(a)') ''

N_world=size(my_lattice%world)
dim_lat=my_lattice%dim_lat
N_cell=product(dim_lat)

allocate(field(1,N_cell))
field=0.0d0


call get_topoplot(field)

Q=get_charge()

N_pat=abs( NINT((Q(1)+Q(2))/pi(4.0d0)) )

write(6,'(a,2x,I3,2x,a)') 'first guess indicates', N_pat, 'patterns found'

order_derivative=1
io_input=open_file_read('input')
call get_parameter(io_input,'input','order_derivative',order_derivative)
call close_file('input',io_input)

allocate(derivative(2*order_derivative+1,N_world,N_cell))

call get_derivative(field,my_lattice,derivative)




allocate(all_patterns(N_pat))
do i=1,N_pat
  call locate_pattern(Nsites,field,N_cell)

  all_patterns(i)%Nsites=Nsites
  allocate(all_patterns(i)%interior(Nsites))
  all_patterns(i)%interior=0

  write(6,'(2(a,2x,I6),a)') 'pattern', i, ' is composed of', Nsites, 'sites'

  call find_interior(all_patterns(i)%Nsites,all_patterns(i)%interior,field,N_cell)

  call find_center(all_patterns(i)%center,all_patterns(i)%Nsites,all_patterns(i)%interior,field,N_cell)

enddo



do i=1,N_pat
  call locate_boundary(Nsites,all_patterns(i)%interior,field,N_cell)
  allocate(all_patterns(i)%boundary(Nsites))
  all_patterns(i)%boundary=0

  call find_boundary(Nsites,all_patterns(i)%boundary,all_patterns(i)%interior,field,derivative,N_cell)

enddo

end subroutine




!
! INTERIOR PART
!



!
! locate the given pattern in the field
!
subroutine locate_pattern(N,field,N_cell)
implicit none
integer, intent(out) :: N
real(kind=8), intent(in) :: field(:,:)
integer, intent(in) :: N_cell
! internal
integer :: i
real(kind=8) :: Q

N=0
do i=1,N_cell
   Q=field(1,i)+field(2,i)
   if (abs(Q).gt.1.0d-8) N=N+1
enddo

end subroutine


!
! locate the finds the interior of the skyrmion
!
subroutine find_interior(Nsites,interior,field,N_cell)
implicit none
integer, intent(in) :: Nsites,N_cell
real(kind=8), intent(in) :: field(:,:)
integer, intent(inout) :: interior(:)
! internal
integer :: i,j
real(kind=8) :: Q

j=0
do i=1,N_cell
   Q=field(1,i)+field(2,i)
   if (abs(Q).gt.1.0d-8) then
     j=j+1
     interior(j)=i
   endif
enddo

end subroutine

!
! locate the center of the pattern interpreted as the barycenter of the field
! x_i= (\sum_j F_j*j)/(\sum_j F_j)
!
subroutine find_center(center,Nsites,interior,field,N_cell)
implicit none
integer, intent(in) :: Nsites,N_cell
real(kind=8), intent(in) :: field(:,:)
integer, intent(in) :: interior(:)
integer, intent(out) :: center
! internal
integer :: i,i_site
real(kind=8) :: Q,sum_Q,sum_Qpos

sum_Q=0.0d0
sum_Qpos=0.0d0
do i=1,Nsites

   i_site=interior(i)
   Q=field(1,i_site)+field(2,i_site)

   sum_Q=sum_Q+Q
   sum_Qpos=sum_Qpos+Q*real(i_site)

enddo

center=NINT(sum_Qpos/sum_Q)

end subroutine


!
! BOUNDARY PART
!

!
! locate the edge of a pattern
!
subroutine locate_boundary(N,interior,field,N_cell)
implicit none
integer, intent(out) :: N
real(kind=8), intent(in) :: field(:,:)
integer, intent(in) :: N_cell,interior(:)
! internal
integer :: i,N_interior,i_site
real(kind=8) :: Qavant,Qapres

N_interior=size(interior)
N=0
do i=1,N_interior
   i_site=interior(i)

   Qavant=field(1,i_site-1)+field(2,i_site-1)
   Qapres=field(1,i_site+1)+field(2,i_site+1)
   if ((abs(Qavant).lt.1.0d-8).or.(abs(Qavant).lt.1.0d-8)) N=N+1
enddo

end subroutine

!
! locate the finds the interior of the skyrmion
!
subroutine find_boundary(Nsites,boundary,interior,field,DF,N_cell)
implicit none
integer, intent(in) :: Nsites,N_cell,interior(:)
real(kind=8), intent(in) :: field(:,:)
type(vec_point), intent(in) :: DF(:,:,:)
integer, intent(inout) :: boundary(:)
! internal
integer :: i,N_interior,i_site,j,shape_DF(3),n_dim,n_vois
real(kind=8) :: Qavant,Qapres,Derivate(2)

N_interior=size(interior)
shape_DF=shape(DF)
n_dim=shape_DF(2)
n_vois=shape_DF(1)
Derivate=0.0d0
j=0
do i=1,N_interior
   i_site=interior(i)

   call calculate_derivative(Derivate,DF(:,:,i_site),n_dim,n_vois)

   Qavant=field(1,i_site-1)+field(2,i_site-1)
   Qapres=field(1,i_site+1)+field(2,i_site+1)
   if ((abs(Qavant).lt.1.0d-8).or.(abs(Qavant).lt.1.0d-8)) then
     j=j+1
     boundary(j)=i_site
   endif
enddo

end subroutine



!
! plot the different patterns
!
subroutine plot_tracking(tag,mode)
use m_convert
use m_io_files_utils
use m_local_energy
implicit none
integer, intent(in) :: tag
type(vec_point), intent(in) :: mode(:)
! internal
integer :: i,N_pat,io_inside,io_border,j,iomp,dim_mode
character(len=50) :: filename_inside,filename_border
real(kind=8) :: Et

filename_inside=convert('pat_inside_',tag,'.dat')
io_inside=open_file_write(filename_inside)

filename_border=convert('pat_border_',tag,'.dat')
io_border=open_file_write(filename_border)

dim_mode=size(mode(1)%w)
N_pat=size(all_patterns)
do i=1,N_pat

  do j=1,all_patterns(i)%Nsites

    iomp=all_patterns(i)%interior(j)
    call local_energy(Et,iomp,mode,dim_mode)

    write(io_inside,'(I8,2x,E20.12E3)') iomp,Et

  enddo

  write(io_inside,'(a)') ''
enddo

dim_mode=size(mode(1)%w)
do i=1,N_pat

  do j=1,size(all_patterns(i)%boundary)

    iomp=all_patterns(i)%boundary(j)
    call local_energy(Et,iomp,mode,dim_mode)

    write(io_border,'(I8,2x,E20.12E3)') iomp,Et

  enddo

  write(io_border,'(a)') ''
enddo

call close_file(filename_inside,io_inside)
call close_file(filename_border,io_border)

end subroutine

end module m_tracker
