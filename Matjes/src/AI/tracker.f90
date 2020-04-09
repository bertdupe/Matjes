module m_tracker
use m_topoplot
use m_topo_commons

! tracking parameters
type pattern
  integer, allocatable :: interior(:)
  integer, allocatable :: boundary(:)
  integer :: center,Nsites
end type

type(pattern), allocatable :: all_patterns(:)

! tracking frequency
integer :: track_freq=1

private
public :: init_tracking, plot_tracking

contains

!
! First guess of what should be tracked
!
subroutine init_tracking(dim_lat)
use m_constants, only : pi
implicit none
integer, intent(in) :: dim_lat(:)
! internal variable
real(kind=8), allocatable :: field(:,:)
integer :: i,N_pat,Nsites,N_cell
real(kind=8) :: Q(5)

N_cell=product(dim_lat)
allocate(field(5,N_cell))
field=0.0d0

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

  write(6,'(2(a,2x,I6),a)') 'pattern', i, 'is composed of', Nsites, 'sites'

  call find_interior(all_patterns(i)%Nsites,all_patterns(i)%interior,field,N_cell)

  call find_center(all_patterns(i)%center,all_patterns(i)%Nsites,all_patterns(i)%interior,field,N_cell)

enddo

do i=1,N_pat
  call locate_boundary(Nsites,all_patterns(i)%interior,field,N_cell)
  allocate(all_patterns(i)%boundary(Nsites))
  all_patterns(i)%boundary=0

  call find_boundary(Nsites,all_patterns(i)%boundary,all_patterns(i)%interior,field,N_cell)

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
subroutine find_boundary(Nsites,boundary,interior,field,N_cell)
implicit none
integer, intent(in) :: Nsites,N_cell,interior(:)
real(kind=8), intent(in) :: field(:,:)
integer, intent(inout) :: boundary(:)
! internal
integer :: i,N_interior,i_site,j
real(kind=8) :: Qavant,Qapres

N_interior=size(interior)
j=0
do i=1,N_interior
   i_site=interior(i)

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
subroutine plot_tracking(tag,mode_E_column,E_line)
use m_convert
use m_io_files_utils
use m_local_energy
use m_derived_types, only : point_shell_Operator,point_shell_mode
implicit none
integer, intent(in) :: tag
type(point_shell_Operator), intent(in) :: E_line(:)
type(point_shell_mode), intent(in) :: mode_E_column(:)
! internal
integer :: i,N_pat,io_inside,io_border,j,iomp
character(len=50) :: filename_inside,filename_border
real(kind=8) :: Et

filename_inside=convert('pat_inside_',tag,'.dat')
io_inside=open_file_write(filename_inside)

filename_border=convert('pat_border_',tag,'.dat')
io_border=open_file_write(filename_border)

N_pat=size(all_patterns)
do i=1,N_pat

  do j=1,all_patterns(i)%Nsites

    iomp=all_patterns(i)%interior(j)
    call local_energy(Et,iomp,mode_E_column(iomp),E_line(iomp))

    write(io_inside,'(I8,2x,E20.12E3)') iomp,Et

  enddo

  write(io_inside,'(a)') ''
enddo

do i=1,N_pat

  do j=1,size(all_patterns(i)%boundary)

    iomp=all_patterns(i)%boundary(j)
    call local_energy(Et,iomp,mode_E_column(iomp),E_line(iomp))

    write(io_border,'(I8,2x,E20.12E3)') iomp,Et

  enddo

  write(io_border,'(a)') ''
enddo

call close_file(filename_inside,io_inside)
call close_file(filename_border,io_border)

end subroutine

end module m_tracker
