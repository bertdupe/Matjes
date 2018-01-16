module m_get_position

interface get_position
module procedure get_position_lattice
end interface get_position

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

subroutine get_position_lattice(position,dim_lat,r,motif)
use m_derived_types
implicit none
real(kind=8), intent(inout) :: position(:,:,:,:,:)
integer, intent(in) :: dim_lat(:)
real(kind=8), intent(in) :: r(:,:)
type (cell), intent(in) :: motif
! internal variables
integer :: i_z,i_y,i_x,i_m,Natom_motif

Natom_motif=size(motif%i_mom)

do i_m=1,Natom_motif
   do i_z=1,dim_lat(3)
      do i_y=1,dim_lat(2)
         do i_x=1,dim_lat(1)

         if (.not.(motif%i_mom(i_m))) cycle

!fix the coordinates. This part takes care of the motif of the structure
         position(1:3,i_x,i_y,i_z,i_m)=r(1,:)*(dble(i_x-1)+motif%pos(i_m,1))+r(2,:)*(dble(i_y-1)+motif%pos(i_m,2))+ &
       r(3,:)*(dble(i_z-1)+motif%pos(i_m,3))

         enddo
      enddo
   enddo
enddo

end subroutine get_position_lattice

end module m_get_position
