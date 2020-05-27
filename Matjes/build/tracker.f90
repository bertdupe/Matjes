module m_tracker
use m_topoplot
use m_topo_commons
use m_derived_types, only : lattice
use m_basic_types, only : vec_point
use m_io_files_utils
use m_get_position

! tracking parameters
type pattern
  integer, allocatable :: interior(:)
  integer, allocatable :: boundary(:)
  integer :: center,Nsites
end type

! derivation of the field of interest
real(kind=8), allocatable :: pos(:,:)

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
use m_io_utils
implicit none
type(lattice), intent(in) :: my_lattice
! internal variable
integer :: i,N_pat,Nsites,N_cell,N_world,dim_lat(3),i_site
! derivative parameter
real(kind=8) :: Q(5),real_vec(3,3)

write(6,'(a)') ''
write(6,'(a)') 'initialization of the tracking procedure'
write(6,'(a)') ''

N_world=size(my_lattice%world)
dim_lat=my_lattice%dim_lat
N_cell=product(dim_lat)
real_vec=transpose(my_lattice%areal)

allocate(field(5,N_cell),pos(3,N_cell))
field=0.0d0
pos=0.0d0

call get_position(pos,'positions.dat ')


call get_topoplot(field)

Q=get_charge()

N_pat=abs( NINT((Q(1)+Q(2))/pi(4.0d0)) )

write(6,'(a,2x,I3,2x,a)') 'first guess indicates', N_pat, 'patterns found'

allocate(all_patterns(N_pat))
do i=1,N_pat

  call locate_pattern(Nsites,field,N_cell)

  all_patterns(i)%Nsites=Nsites
  allocate(all_patterns(i)%interior(Nsites))
  all_patterns(i)%interior=0

  write(6,'(2(a,2x,I6),a)') 'pattern  ', i, '  is composed of  ', Nsites, '  sites'

  call find_interior(all_patterns(i)%Nsites,all_patterns(i)%interior,field,N_cell)

  call find_center(all_patterns(i)%center,all_patterns(i)%Nsites,all_patterns(i)%interior,field,N_cell)

  write(6,'(2(a,2x,I6))') 'center of ', i, ' at position ', all_patterns(i)%center

enddo

do i=1,N_pat

  call find_boundary(Nsites,i,field,N_cell,real_vec)

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
real(kind=8) :: Q,sum_Q,sum_Qpos(3),test

sum_Q=0.0d0
sum_Qpos=0.0d0
do i=1,Nsites

   i_site=interior(i)
   Q=field(1,i_site)+field(2,i_site)

   sum_Q=sum_Q+Q
   sum_Qpos=sum_Qpos+Q*pos(:,i_site)

enddo

test=sum((pos(:,1)-sum_Qpos/sum_Q)**2)
center=1
do i=1,Nsites
  i_site=interior(i)
  if (sum((pos(:,i_site)-sum_Qpos/sum_Q)**2).lt.test) then
    test=sum((pos(:,i_site)-sum_Qpos/sum_Q)**2)
    center=i_site
  endif
enddo

end subroutine

!
! locate the finds the interior of the skyrmion
!
subroutine find_boundary(Nsites,i_pat,field,N_cell,real_vec)
use m_vector, only : norm
use m_matrix, only : invert
use m_sort
use m_constants, only : pi
use m_envelope
implicit none
integer, intent(in) :: Nsites,N_cell,i_pat
real(kind=8), intent(in) :: field(:,:),real_vec(:,:)
! internal
integer :: i,i_site,i_center
! test of the boundary
integer, allocatable :: ind(:),vertex(:)
real(kind=8), allocatable :: x(:),y(:)
integer :: N_interior,nvert
! position of the center
real(kind=8) :: P0(3),vec(3)

N_interior=size(all_patterns(i_pat)%interior)
! find the invert of the areal matrix to pass from the
allocate(x(N_interior),y(N_interior))
x=0.0d0
y=0.0d0
i_center=all_patterns(i_pat)%center
P0=pos(:,i_center)

do i=1,N_interior
   i_site=all_patterns(i_pat)%interior(i)

   vec=pos(:,i_site)-P0
   x(i)=vec(1)
   y(i)=vec(2)

enddo

allocate(ind(N_interior),vertex(N_interior))
ind=0
vertex=0

call envelope(x, y, N_interior, vertex, nvert, ind)

allocate(all_patterns(i_pat)%boundary(nvert))
do i=1,nvert
  i_site=vertex(i)
  all_patterns(i_pat)%boundary(i)=all_patterns(i_pat)%interior(i_site)
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

    write(io_inside,'(4(2x,E20.12E3))') pos(:,iomp),Et

  enddo

  write(io_inside,'(a)') ''
enddo

dim_mode=size(mode(1)%w)
do i=1,N_pat

  do j=1,size(all_patterns(i)%boundary)

    iomp=all_patterns(i)%boundary(j)

    call local_energy(Et,iomp,mode,dim_mode)

    write(io_border,'(4(2x,E20.12E3))') pos(:,iomp),Et

  enddo

  write(io_border,'(a)') ''
enddo

call close_file(filename_inside,io_inside)
call close_file(filename_border,io_border)

end subroutine

end module m_tracker
