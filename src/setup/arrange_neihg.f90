module m_arrange_neigh
use m_rotation, only : rotation_axis,rotate
use m_vector, only : calc_ang,cross
! get the symmetry operations used to rotate the neighbors
! the reference direction is the x direction
type rotation_operation
  real(kind=8) :: rotation_axis(3)
  real(kind=8) :: rotation_angle
end type rotation_operation

type shell_maps
  type(rotation_operation), allocatable, dimension(:) :: bond
end type shell_maps

type(shell_maps), allocatable, dimension(:) :: map


!interface arrange_neigh
!  module procedure arrange_neigh_DM
!end interface arrange_neigh

private
public :: arrange_DM,get_map,fill_map,get_rotation

contains

!
! define the map matrix
!
subroutine get_map(indexNN)
implicit none
integer, intent(in) :: indexNN(:,:)
! internal
integer :: i,j

if(allocated(map))then
    WRITE(*,*) "WARNING, GET_MAP more than once CALLED"    
    !I only added this in order to hack the tight-binding neighbors
    deallocate(map) !added ugly
endif
allocate( map(size(indexNN,1)) )
do i=1,size(indexNN,1)
  allocate( map(i)%bond( indexNN(i,1) ) )
  do j=1,indexNN(i,1)
    map(i)%bond(j)%rotation_axis=0.0d0
    map(i)%bond(j)%rotation_angle=0.0d0
  enddo
enddo

end subroutine get_map

!
! fill the map matrix
!
subroutine fill_map(i_shell,i_bond,vec)
use m_rotation, only : rotation_axis
use m_vector, only : calc_ang,norm
implicit none
integer, intent(in) :: i_shell, i_bond
real(kind=8), intent(in) :: vec(:)
! internal
real(kind=8) :: axis(3)

axis=0.0d0

map(i_shell)%bond(i_bond)%rotation_angle=calc_ang((/1.0d0,0.0d0,0.0d0/),vec)
axis=rotation_axis((/1.0d0,0.0d0,0.0d0/),vec)
if (norm(axis).lt.1.0d-8) then
  !
  ! the bound is along x or -x depending on the angle
  !
  map(i_shell)%bond(i_bond)%rotation_axis=(/0.0d0,0.0d0,1.0d0/)
else
  map(i_shell)%bond(i_bond)%rotation_axis=axis
endif

end subroutine fill_map


!
! arrange the DMI vector based on the rotation of the bounds
!

subroutine arrange_DM(DM_vector,n_DMI)
use m_constants, only : pi
implicit none
real(kind=8) , intent(inout) :: DM_vector(:,:,:)
integer, intent(in) :: n_DMI
! internal
real(kind=8), allocatable, dimension(:,:,:) :: int_DM_vector
real(kind=8) :: cos_angle_bond, axe_rotation_local(3),bound_direction(3)
integer :: shape_DM(3),avant
integer :: i_dm,i_bond,n_bound,i_test,counter

shape_DM=shape(DM_vector)
allocate(int_DM_vector(shape_DM(1),shape_DM(2),shape_DM(3)))
int_DM_vector=0.0d0

avant=0
counter=0
do i_dm=1,n_DMI
  n_bound=size(map(i_dm)%bond)
  do i_bond=1,n_bound

    bound_direction=0.0d0
    ! map(i_dm)%bond(i_bond)%rotation_axis is the rotation axis and map(i_dm)%bond(i_bond)%rotation_angle is the rotation angle taht transforms x direction
    ! to the nearst neighbor bound n°i_bound
    ! bound_direction is the direction of the neighbor i_bond
    call rotate((/1.0d0,0.0d0,0.0d0/),map(i_dm)%bond(i_bond)%rotation_axis,map(i_dm)%bond(i_bond)%rotation_angle,bound_direction)

    do i_test=1,n_bound

      ! the bound must be perpendicular to the DM vector
      cos_angle_bond=dot_product(DM_vector(i_test+avant,:,1),bound_direction)
      ! the rotation between the bound and the DM vector muct be counter-clockwise
      !  axe_rotation_local must be along the z vector
      axe_rotation_local=rotation_axis(bound_direction,DM_vector(i_test+avant,:,1))

!
! that will not work with out-of-plane DMI
!
      if ((abs(cos_angle_bond).lt.1.0d-8).and.(dot_product((/0.0d0,0.0d0,1.0d0/),axe_rotation_local).gt.0.0d0)) then
        counter=counter+1
        int_DM_vector(counter,:,1)=DM_vector(i_test+avant,:,1)
        exit
      endif

    enddo

  enddo
  avant=avant+n_bound
enddo

DM_vector=int_DM_vector

write(6,'(a)') ''
write(6,'(a)') 'reorganization of the DM vectors'
write(6,'(a)') '1 DM must not be perpendicular to x direction (ideally along x)'
do i_dm=1,size(DM_vector,1)
   write(6,'(3(f8.5,2x))') DM_vector(i_dm,:,1)
enddo
write(6,'(a)') ''

end subroutine arrange_DM



!
! get the rotation matrix to transform the operator along the x direction
! to an arbitrary bond direction
! we calculate the image of the unit vectors via the rotation of -angle since we want the matrix that transform (x,y,z) into (x,y,z)'
!

subroutine get_rotation(i_shell,i_bond,h_dum)
use m_rotation, only : rotate
implicit none
integer, intent(in) :: i_shell,i_bond
real(kind=8) , intent(out) :: h_dum(:,:)
! internal

call rotate((/1.0d0,0.0d0,0.0d0/),map(i_shell)%bond(i_bond)%rotation_axis,map(i_shell)%bond(i_bond)%rotation_angle,h_dum(:,1))
call rotate((/0.0d0,1.0d0,0.0d0/),map(i_shell)%bond(i_bond)%rotation_axis,map(i_shell)%bond(i_bond)%rotation_angle,h_dum(:,2))
call rotate((/0.0d0,0.0d0,1.0d0/),map(i_shell)%bond(i_bond)%rotation_axis,map(i_shell)%bond(i_bond)%rotation_angle,h_dum(:,3))

end subroutine get_rotation

end module m_arrange_neigh
