module m_operator_pointer_utils
use m_derived_types, only : operator_real

interface associate_pointer
   module procedure associate_pointer_real_3D
end interface associate_pointer

private
public :: associate_pointer,dissociate
contains


subroutine associate_pointer_real_3D(point,static_target,my_lattice,tableNN,indexNN)
use m_derived_types, only : lattice,operator_real,shell_Ham
use m_get_position
implicit none
type(operator_real), intent(inout) :: point
!type(Coeff_Ham), intent(in) :: static_target
type(shell_Ham), target, intent(in) :: static_target(:)
type(lattice), intent(in) :: my_lattice
integer, intent(in) :: tableNN(:,:,:,:,:,:) !!tableNN(4,N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),count(my_motif%i_mom)
integer, intent(in) :: indexNN(:)
! internal
integer :: i_x,i_y,i_z,i_m,ipos_1,ipos_2,v_x,v_y,v_z,v_m,i_voisin,i_shell
integer :: Nspin,all_size(4),shape_tableNN(6),avant,Ilat(4)
integer :: N_shell
real(kind=8) :: test1,test2

all_size=shape(my_lattice%l_modes)
Nspin=product(all_size)
shape_tableNN=shape(tableNN)
N_shell=size(static_target)

do i_m=1,shape_tableNN(6)
  do i_z=1,shape_tableNN(5)
    do i_y=1,shape_tableNN(4)
      do i_x=1,shape_tableNN(3)

        Ilat=(/i_x,i_y,i_z,i_m/)
        ipos_2=get_position_ND_to_1D(Ilat,all_size)

        ! anisotropy

        point%value(ipos_2,ipos_2)%Op_loc=>static_target(1)%atom(1)%H(:,:)

        avant=0
        do i_shell=1,N_shell-1

          test1=sum(abs(static_target(i_shell+1)%atom(1)%H(:,:)))

          if (test1.lt.1.0d-8) cycle

          do i_voisin=1,indexNN(i_shell)

          v_x=tableNN(1,i_voisin+avant,i_x,i_y,i_z,i_m)
          v_y=tableNN(2,i_voisin+avant,i_x,i_y,i_z,i_m)
          v_z=tableNN(3,i_voisin+avant,i_x,i_y,i_z,i_m)
          v_m=tableNN(4,i_voisin+avant,i_x,i_y,i_z,i_m)

          Ilat=(/v_x,v_y,v_z,v_m/)
          ipos_1=get_position_ND_to_1D(Ilat,all_size)

          test2=sum(abs(static_target(i_shell+1)%atom(i_voisin)%H(:,:)))

          if (test2.gt.1.0d-8) then
             point%value(ipos_1,ipos_2)%Op_loc=>static_target(i_shell+1)%atom(i_voisin)%H(:,:)
          else
             point%value(ipos_1,ipos_2)%Op_loc=>static_target(i_shell+1)%atom(1)%H(:,:)
          endif

!#ifdef CPP_DEBUG
!          write(*,*) v_x,v_y,v_z,v_m
!          write(*,*) i_x,i_y,i_z,i_m
!          write(*,*) static_target(i_shell+1)%atom(i_voisin)%H(1,:)
!          write(*,*) static_target(i_shell+1)%atom(i_voisin)%H(2,:)
!          write(*,*) static_target(i_shell+1)%atom(i_voisin)%H(3,:)
!          pause
!#endif

          enddo
          avant=avant+indexNN(i_shell)

        enddo

      enddo
    enddo
  enddo
enddo

end subroutine associate_pointer_real_3D

subroutine dissociate(matrix,N,M)
implicit none
integer, intent(in) :: N,M
type(operator_real), intent(inout) :: matrix
! internal
integer :: i,j

do i=1,N
   do j=1,M

   nullify(matrix%value(i,j)%Op_loc)

   enddo
enddo

end subroutine dissociate

end module m_operator_pointer_utils
