module m_get_position

  interface get_position
      module procedure get_position_lattice
  end interface get_position

private
public :: get_position,get_position_ND_to_1D
contains

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
use m_derived_types, only : cell
implicit none
real(kind=8), intent(inout) :: pos(:,:,:,:,:)
integer, intent(in) :: dim_lat(:)
real(kind=8), intent(in) :: r(:,:)
type (cell), intent(in) :: motif
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
use m_derived_types, only : cell
implicit none
type(cell), intent(in) :: motif
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
