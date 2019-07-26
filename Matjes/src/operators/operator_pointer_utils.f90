module m_operator_pointer_utils

interface associate_pointer
   module procedure A_vecpoint1D_vecdimn,A_Opreal_shellHam1D,A_Opreal_real2D,associate_line_target
end interface associate_pointer

interface dissociate
   module procedure dissociate_OpReal_2D,diss_point_shell_Operator_1D,diss_point_shell_mode_1D,dissociate_basicOpReal_2D
end interface dissociate

private
public :: associate_pointer,dissociate,associate_line_target

contains
!
! routine that associates a (:) pointer to a 1D matrix
! this is typically used when the electromagnetic field is associated. At least the amplitude
!
subroutine A_vecpoint1D_vecdimn(point,static_target,tableNN)
use m_derived_types, only : vec_point,vec_dim_n
use m_get_position
implicit none
integer, intent(in) :: tableNN(:,:,:,:,:,:) !!tableNN(4,N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),count(my_motif%i_mom)
type(vec_point), intent(inout) :: point(:)
!type(Coeff_Ham), intent(in) :: static_target
type(vec_dim_n), target, intent(in) :: static_target
! internal
integer :: i,Nspin,N_size_target,i_x,i_y,i_z,i_m
integer :: ilat(4),shape_tableNN(6)

Nspin=size(point)
N_size_target=size(static_target%w)
shape_tableNN=shape(tableNN)

do i_m=1,shape_tableNN(6)
  do i_z=1,shape_tableNN(5)
    do i_y=1,shape_tableNN(4)
      do i_x=1,shape_tableNN(3)

         Ilat=(/i_x,i_y,i_z,i_m/)
         i=get_position_ND_to_1D(Ilat,shape_tableNN(3:6))

         point(i)%w=>static_target%w

      enddo
    enddo
  enddo
enddo

end subroutine A_vecpoint1D_vecdimn

!
! routine that associates a (:,:) pointer to a Hamiltonian matrix
! typical assocication for the Hamiltonian
!
subroutine A_Opreal_shellHam1D(point,static_target,my_lattice,tableNN,indexNN)
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
integer :: Nspin,all_size(4),shape_tableNN(6),avant,Ilat(4),line_index
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

        line_index=1

        ! anisotropy

        point%value(line_index,ipos_2)%Op_loc=>static_target(1)%atom(1)%H(:,:)
        point%line(line_index,ipos_2)=ipos_2

        avant=0
        do i_shell=1,N_shell-1

          do i_voisin=1,indexNN(i_shell)

          line_index=line_index+1

          v_x=tableNN(1,i_voisin+avant,i_x,i_y,i_z,i_m)
          v_y=tableNN(2,i_voisin+avant,i_x,i_y,i_z,i_m)
          v_z=tableNN(3,i_voisin+avant,i_x,i_y,i_z,i_m)
          v_m=tableNN(4,i_voisin+avant,i_x,i_y,i_z,i_m)

          Ilat=(/v_x,v_y,v_z,v_m/)
          ipos_1=get_position_ND_to_1D(Ilat,all_size)

          point%line(line_index,ipos_2)=ipos_1

! if the matrix is symmetric (exchnage only) then there is no point is storing
! as much energy matrices as in the number of atom in the shell.
!
! Therefore if the energy for the atom i_voisin in shell i_shell is 0 it does not mean that Hamiltonian of this shell
! should be 0. One has to check the Hamiltonian for the 1st atom in the shell
          test2=sum(abs(static_target(i_shell+1)%atom(i_voisin)%H(:,:)))

          if (test2.gt.1.0d-8) then
             point%value(line_index,ipos_2)%Op_loc=>static_target(i_shell+1)%atom(i_voisin)%H(:,:)
          else
             point%value(line_index,ipos_2)%Op_loc=>static_target(i_shell+1)%atom(1)%H(:,:)
          endif

!#ifdef CPP_DEBUG
!          write(*,*) v_x,v_y,v_z,v_m
!          write(*,*) i_x,i_y,i_z,i_m
!          write(*,*) ipos_1,ipos_2,line_index
!          write(*,*) static_target(i_shell+1)%atom(i_voisin)%H(1,:)
!          write(*,*) static_target(i_shell+1)%atom(i_voisin)%H(2,:)
!          write(*,*) static_target(i_shell+1)%atom(i_voisin)%H(3,:)
!          write(*,*) point%value(line_index,ipos_2)%Op_loc(1,:)
!          write(*,*) point%value(line_index,ipos_2)%Op_loc(2,:)
!          write(*,*) point%value(line_index,ipos_2)%Op_loc(3,:)
!          pause
!#endif

          enddo
          avant=avant+indexNN(i_shell)

        enddo

      enddo
    enddo
  enddo
enddo

end subroutine A_Opreal_shellHam1D

!
! routine that associates a (:,:) pointer to a 2D matrix of reals
! typical assocication for the shell decomposition of the Hamiltonian
!
subroutine A_Opreal_real2D(point,static_target,my_lattice,tableNN,Nvoisin,avant)
use m_derived_types, only : lattice,operator_real,shell_Ham
use m_get_position
implicit none
type(operator_real), intent(inout) :: point
real(kind=8), target, intent(in) :: static_target(:,:)
type(lattice), intent(in) :: my_lattice
integer, intent(in) :: tableNN(:,:,:,:,:,:) !!tableNN(4,N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),count(my_motif%i_mom)
integer, intent(in) :: Nvoisin
integer, intent(in), optional :: avant
! internal
integer :: i_x,i_y,i_z,i_m,ipos_1,ipos_2,v_x,v_y,v_z,v_m,i_voisin,i_shell
integer :: line_index,advanced
integer :: Nspin,all_size(4),shape_tableNN(6),Ilat(4)
integer :: N_shell
real(kind=8) :: test1,test2

all_size=shape(my_lattice%l_modes)
Nspin=product(all_size)
shape_tableNN=shape(tableNN)
N_shell=size(static_target)

if (.not.present(avant)) then
   advanced=0
else
   advanced=avant
endif

do i_m=1,shape_tableNN(6)
  do i_z=1,shape_tableNN(5)
    do i_y=1,shape_tableNN(4)
      do i_x=1,shape_tableNN(3)

        Ilat=(/i_x,i_y,i_z,i_m/)
        ipos_2=get_position_ND_to_1D(Ilat,all_size)

        line_index=1

        if (advanced.eq.0) then
        ! anisotropy

           point%value(line_index,ipos_2)%Op_loc=>static_target

           point%line(line_index,ipos_2)=ipos_2

        else

          do i_voisin=1,Nvoisin

            line_index=line_index+1

            v_x=tableNN(1,i_voisin+advanced,i_x,i_y,i_z,i_m)
            v_y=tableNN(2,i_voisin+advanced,i_x,i_y,i_z,i_m)
            v_z=tableNN(3,i_voisin+advanced,i_x,i_y,i_z,i_m)
            v_m=tableNN(4,i_voisin+advanced,i_x,i_y,i_z,i_m)

            Ilat=(/v_x,v_y,v_z,v_m/)
            ipos_1=get_position_ND_to_1D(Ilat,all_size)

            point%value(line_index,ipos_2)%Op_loc=>static_target
            point%line(line_index,ipos_2)=ipos_1

            write(*,*) point%value(line_index,ipos_2)%Op_loc

          enddo
        endif

      enddo
    enddo
  enddo
enddo

end subroutine A_Opreal_real2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! optimize the calculation time
! associate each non-zero elements in the columns of the operator to the good lines in the local modes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine associate_line_target(line,mode,column,target)
use m_derived_types, only : operator_real,vec_point,point_shell_Operator,point_shell_mode
implicit none
type(vec_point),target,intent(in) :: mode(:)
type(operator_real),target,intent(in) :: target
type(point_shell_Operator), intent(inout) :: column(:)
type(point_shell_mode), intent(inout) :: line(:)
!internal energy
integer :: i,j,shape_target(2),size_mode,pos_cell

shape_target=shape(target%value)
size_mode=size(mode)

if (shape_target(2).ne.size_mode) stop 'error 1 in associate_line_target'

do i=1,size_mode

   do j=1,shape_target(1)

      pos_cell=target%line(j,i)

      line(i)%shell(j)%w=>mode(pos_cell)%w
      column(i)%shell(j)%Op_loc=>target%value(j,i)%Op_loc

   enddo

enddo

end subroutine associate_line_target

!!!!!!!!!!!!!!!!!
!
! Dissociate the pointers from anything
!
!!!!!!!!!!!!!!!!!
subroutine dissociate_basicOpReal_2D(matrix,N,M)
use m_basic_types, only : Op_real
implicit none
integer, intent(in) :: N,M
type(Op_real), intent(out) :: matrix(N,M)
! internal
integer :: i,j

do j=1,M
   do i=1,N

   nullify(matrix(i,j)%Op_loc)

   enddo
enddo

end subroutine dissociate_basicOpReal_2D

!!!!!!!!!!!!!!!!!
!
! Dissociate the pointers from anything
!
!!!!!!!!!!!!!!!!!
subroutine dissociate_OpReal_2D(matrix)
use m_derived_types, only : operator_real
implicit none
type(operator_real), intent(inout) :: matrix
! internal
integer :: i,j,size_matrix(2)

size_matrix=size(matrix%value)

do j=1,size_matrix(2)
   do i=1,size_matrix(1)

   nullify(matrix%value(i,j)%Op_loc)

   enddo
enddo

end subroutine dissociate_OpReal_2D

!!!!!!!!!!!!!!!!!
!
! Dissociate the line of pointers from anything
!
!!!!!!!!!!!!!!!!!
subroutine diss_point_shell_Operator_1D(matrix,M,N)
use m_derived_types, only : point_shell_Operator
implicit none
type(point_shell_Operator), intent(inout) :: matrix(:)
! internal
integer :: i,j,N,M

do j=1,N
  do i=1,M
    nullify(matrix(j)%shell(i)%Op_loc)
  enddo

enddo

end subroutine diss_point_shell_Operator_1D

!!!!!!!!!!!!!!!!!
!
! Dissociate the line of pointers from anything
!
!!!!!!!!!!!!!!!!!
subroutine diss_point_shell_mode_1D(matrix,M,N)
use m_derived_types, only : point_shell_mode
implicit none
type(point_shell_mode), intent(inout) :: matrix(:)
! internal
integer :: i,j,N,M

do j=1,N
  do i=1,M
    nullify(matrix(j)%shell(i)%w)
  enddo

enddo

end subroutine diss_point_shell_mode_1D

end module m_operator_pointer_utils
