module m_get_position
use m_io_files_utils

  interface get_position
      module procedure get_position_lattice
      module procedure get_position_file_4d
      module procedure get_position_file_1d
  end interface get_position

private
public :: get_position,get_position_ND_to_1D,calculate_distances
contains

!
! calculate the shorter between 2 points on the lattice (the matrix of the r to calculate the 1/r)
!

subroutine calculate_distances(dist,pos,r,dim_lat,boundary)
use m_vector, only : norm
implicit none
real(kind=8), intent(inout) :: dist(:,:)
real(kind=8), intent(in) :: pos(:,:),r(:,:)
integer, intent(in) :: dim_lat(:)
logical, intent(in) :: boundary(:)
! internal
integer :: i,k,l,m,N
real(kind=8) :: shorter_distance(3),test(3),N1,N2,N3

N1=real(dim_lat(1))
N2=real(dim_lat(2))
N3=real(dim_lat(3))
N=product(dim_lat)
#ifdef CPP_OPENMP
!$OMP parallel do private(shorter_distance,test) default(shared)
#endif

do i=1,N

  test=pos(:,i)

  do k=-1,1
    do l=-1,1
      do m=-1,1

!    direction 1
        if (boundary(1)) shorter_distance=test+real(k)*r(:,1)*N1

!     direction 2
        if (boundary(2)) shorter_distance=shorter_distance+real(l)*r(:,2)*N2

!      direction 3
        if (boundary(3)) shorter_distance=shorter_distance+real(m)*r(:,3)*N3

        if (norm(shorter_distance).lt.norm(test)) test=shorter_distance

      enddo
    enddo
  enddo

  dist(:,i)=test

enddo
#ifdef CPP_OPENMP
!$OMP end parallel do
#endif
end subroutine calculate_distances




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine that reads the position from file
! input:
! dim_lat(3): dimension of the lattice in the x, y and z directions
!
! outpout:
! position: position of the sites in cartesian coordinates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_position_file_4d(pos,fname)
implicit none
real(kind=8), intent(inout) :: pos(:,:,:,:,:)
character(len=*), intent(in) :: fname
! internal variables
integer :: dim_lat(5),io
integer :: i_z,i_y,i_x,i_m,j_lat

dim_lat= shape(pos)

! Read the configurations
io=open_file_read(fname)

do i_m=1,dim_lat(5)
   do i_z=1,dim_lat(4)
      do i_y=1,dim_lat(3)
         do i_x=1,dim_lat(2)
           read(io,*) (pos(j_lat,i_x,i_y,i_z,i_m),j_lat=1,3)
         enddo
      enddo
   enddo
enddo

call close_file(fname,io)

end subroutine get_position_file_4d

subroutine get_position_file_1d(pos,fname)
implicit none
real(kind=8), intent(inout) :: pos(:,:)
character(len=*), intent(in) :: fname
! internal variables
integer :: i,j_lat,io

! Read the configurations
io=open_file_read(fname)

do i=1,size(pos,2)
   read(io,*) (pos(j_lat,i),j_lat=1,3)
enddo

call close_file(fname,io)

end subroutine get_position_file_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine that gives out the position of the of the sites
! input:
! dim_lat(3): dimension of the lattice in the x, y and z directions
! r(3,3): lattice vectors in cartesian
! my_motif: position of the atoms in the motif
!
! outpout:
! position: position of the sites in cartesian coordinates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_position_lattice(pos,dim_lat,r,motif)
use m_derived_types, only : t_cell
implicit none
real(kind=8), intent(inout) :: pos(:,:,:,:,:)
integer, intent(in) :: dim_lat(:)
real(kind=8), intent(in) :: r(:,:)
type(t_cell), intent(in) :: motif
! internal variables
integer :: i_z,i_y,i_x,i_m,Natom_motif

Natom_motif=size(motif%atomic)

do i_m=1,Natom_motif
   do i_z=1,dim_lat(3)
      do i_y=1,dim_lat(2)
         do i_x=1,dim_lat(1)

         if (motif%atomic(i_m)%moment.lt.1.0d-8) cycle

!fix the coordinates. This part takes care of the motif of the structure
         pos(1:3,i_x,i_y,i_z,i_m)=get_internal_to_cart(r,motif,i_x,i_y,i_z,i_m)

         enddo
      enddo
   enddo
enddo

end subroutine get_position_lattice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function that returns that maps a ND lattice (i_x,i_y,i_z,i_m...) on a 1D lattice
! input:
! Ilat(:): position on the ND lattice
! N(:): dimension of the ND lattice
!
! outpout:
! get_position_4D_to_1D: position in the 1D lattice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive function get_position_ND_to_1D(Ilat,N) result(k)
use m_derived_types
implicit none
integer, intent(in) :: N(:),Ilat(:)
integer :: k
! internal variables
integer :: Num,Ncoord

Num=size(N)
Ncoord=size(ilat)

if (Num.eq.Ncoord) k=0

if (Num.eq.1) then
   k=ilat(Num)
else
   k=(ilat(Num)-1)*product(N(1:Num-1))+get_position_ND_to_1D(Ilat,N(1:Num-1))
endif

end function get_position_ND_to_1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subfunction that gives cartesian corrdinates as a function of the internal coordinates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function get_internal_to_cart(r,motif,i_x,i_y,i_z,i_m)
use m_derived_types, only : t_cell
implicit none
type(t_cell), intent(in) :: motif
real(kind=8), intent(in) :: r(:,:)
integer, intent(in) :: i_x,i_y,i_z,i_m
! output
real(kind=8) :: get_internal_to_cart(3)
! internal

get_internal_to_cart=0.0d0

get_internal_to_cart=r(1,:)*(dble(i_x-1)+motif%atomic(i_m)%position(1))+r(2,:)*(dble(i_y-1)+motif%atomic(i_m)%position(2))+ &
       r(3,:)*(dble(i_z-1)+motif%atomic(i_m)%position(3))

end function

end module m_get_position
