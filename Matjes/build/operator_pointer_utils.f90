module m_operator_pointer_utils

interface associate_pointer
   module procedure A_vecpoint1D_vecdimn,A_Opreal_shellHam1D,A_Opreal_real2D,associate_line_target
   module procedure associate_mode_name,associate_mode_real2D_name
end interface associate_pointer

interface dissociate
   module procedure dissociate_OpReal_2D,diss_point_shell_Operator_1D,diss_point_shell_mode_1D,dissociate_basicOpReal_2D
   module procedure diss_point_shell_1D_Operator_1D,diss_point_shell_1D_mode_1D,diss_vec_point_1D
end interface dissociate

private
public :: associate_pointer,dissociate,associate_line_target

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! routine that associates a (:) pointer to part of the first dimension of a 2D matrix
! this is typically used when you only want to get the magnetic moments or the integrator fields
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine associate_mode_real2D_name(point,static_target,name,i_name)
use m_derived_types, only : vec_point
use m_lattice, only : my_order_parameters
implicit none
type(vec_point), intent(inout) :: point(:)
logical, intent(inout) :: i_name
real(kind=8), target, intent(in) :: static_target(:,:)
character(len=*), intent(in) :: name
! internal
integer :: i,size_point,size_target,N_orders,N_order_found,istart,iend
integer, allocatable :: position_found(:,:)

size_point=size(point)
size_target=size(static_target,2)

if (size_point.ne.size_target) stop 'error in associate_mode_real2D_name - DIM unequal'

! function at the end of the module
N_order_found=numer_order_param(name)
if (N_order_found.ne.0) then
   i_name=.true.
else
   return
endif

allocate(position_found(3,N_order_found))
position_found=0

! function at the end of the module
position_found=find_position_order(N_order_found,name)

istart=position_found(2,1)
iend=position_found(3,N_order_found)

do i=1,size_point
   point(i)%w=>static_target(istart:iend,i)
enddo

end subroutine associate_mode_real2D_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! routine that associates a (:) pointer to a part of another pointer
! this is typically used when you only want to get the magnetic moments or the local modes and not all modes
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine associate_mode_name(point,static_target,name,i_name)
use m_derived_types, only : vec_point
use m_lattice, only : my_order_parameters
implicit none
type(vec_point), intent(inout) :: point(:)
logical, intent(inout) :: i_name
type(vec_point), intent(in) :: static_target(:)
character(len=*), intent(in) :: name
! internal
integer :: i,size_point,size_target,N_orders,N_order_found,istart,iend
integer, allocatable :: position_found(:,:)

size_point=size(point)
size_target=size(static_target)

if (size_point.ne.size_target) stop 'error in associate_mode_name - DIM unequal'

! function at the end of the module
N_order_found=numer_order_param(name)
if (N_order_found.ne.0) then
   i_name=.true.
else
   return
endif

allocate(position_found(3,N_order_found))
position_found=0

! function at the end of the module
position_found=find_position_order(N_order_found,name)

istart=position_found(2,1)
iend=position_found(3,N_order_found)

do i=1,size_point
   point(i)%w=>static_target(i)%w(istart:iend)
enddo

end subroutine associate_mode_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! routine that associates a (:) pointer to a 1D matrix
! this is typically used when the electromagnetic field is associated. At least the amplitude
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! routine that associates a (:,:) pointer to a Hamiltonian matrix
! typical assocication for the Hamiltonian
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

        point%value(line_index,ipos_2)%Op_loc=>static_target(1)%atom(1)%H
        point%line(line_index,ipos_2)=ipos_2

        avant=0
        do i_shell=1,N_shell-1

          do i_voisin=1,size(static_target(i_shell+1)%atom)

          line_index=line_index+1

          v_x=tableNN(1,i_voisin+avant,i_x,i_y,i_z,i_m)
          v_y=tableNN(2,i_voisin+avant,i_x,i_y,i_z,i_m)
          v_z=tableNN(3,i_voisin+avant,i_x,i_y,i_z,i_m)
          v_m=tableNN(4,i_voisin+avant,i_x,i_y,i_z,i_m)

          Ilat=(/v_x,v_y,v_z,v_m/)
          ipos_1=get_position_ND_to_1D(Ilat,all_size)

          point%line(line_index,ipos_2)=ipos_1

          point%value(line_index,ipos_2)%Op_loc=>static_target(i_shell+1)%atom(i_voisin)%H

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
          avant=avant+size(static_target(i_shell+1)%atom)

        enddo

      enddo
    enddo
  enddo
enddo

end subroutine A_Opreal_shellHam1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! routine that associates a (:,:) pointer to a 2D matrix of reals
! typical assocication for the shell decomposition of the Hamiltonian
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine A_Opreal_real2D(point,static_target,my_lattice,tableNN,Nvoisin,avant)
use m_derived_types, only : lattice,operator_real,shell_Ham
use m_get_position
implicit none
type(operator_real), intent(inout) :: point
real(kind=8), target, intent(in) :: static_target(:,:)
type(lattice), intent(in) :: my_lattice
integer, intent(in) :: tableNN(:,:,:,:,:,:) !!tableNN(4,N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),count(my_motif%i_mom)
integer, intent(in) :: Nvoisin
integer, intent(in) :: avant
! internal
integer :: i_x,i_y,i_z,i_m,ipos_1,ipos_2,v_x,v_y,v_z,v_m,i_voisin
integer :: Nspin,all_size(4),shape_tableNN(6),Ilat(4)
integer :: N_shell

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


          do i_voisin=1,Nvoisin

            v_x=tableNN(1,i_voisin+avant,i_x,i_y,i_z,i_m)
            v_y=tableNN(2,i_voisin+avant,i_x,i_y,i_z,i_m)
            v_z=tableNN(3,i_voisin+avant,i_x,i_y,i_z,i_m)
            v_m=tableNN(4,i_voisin+avant,i_x,i_y,i_z,i_m)

            Ilat=(/v_x,v_y,v_z,v_m/)
            ipos_1=get_position_ND_to_1D(Ilat,all_size)

            point%value(i_voisin,ipos_2)%Op_loc=>static_target
            point%line(i_voisin,ipos_2)=ipos_1

          enddo

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
subroutine dissociate_OpReal_2D(matrix,M,N)
use m_derived_types, only : operator_real
implicit none
type(operator_real), intent(inout) :: matrix
! internal
integer :: i,j,M,N

do j=1,M
   do i=1,N

   nullify(matrix%value(N,M)%Op_loc)

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
integer, intent(in) :: N,M
! internal
integer :: i,j

do j=1,N
  do i=1,M
    nullify(matrix(j)%shell(i)%Op_loc)
  enddo

enddo

end subroutine diss_point_shell_Operator_1D

!!!!!!!!!!!!!!!!!
!
! Dissociate the line of pointers when the dimensions are not aloways the same from anything
!
!!!!!!!!!!!!!!!!!
subroutine diss_point_shell_1D_Operator_1D(matrix,M,N)
use m_derived_types, only : point_shell_Operator
implicit none
type(point_shell_Operator), intent(inout) :: matrix(:)
integer, intent(in) :: M(:),N
! internal
integer :: i,j

do j=1,N
  do i=1,M(j)
    nullify(matrix(j)%shell(i)%Op_loc)
  enddo

enddo

end subroutine diss_point_shell_1D_Operator_1D

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

!!!!!!!!!!!!!!!!!
!
! Dissociate the line of pointers with different dimension from anything
!
!!!!!!!!!!!!!!!!!
subroutine diss_point_shell_1D_mode_1D(matrix,M,N)
use m_derived_types, only : point_shell_mode
implicit none
type(point_shell_mode), intent(inout) :: matrix(:)
integer, intent(in) :: M(:),N
! internal
integer :: i,j

do j=1,N
  do i=1,M(j)
    nullify(matrix(j)%shell(i)%w)
  enddo

enddo

end subroutine diss_point_shell_1D_mode_1D

!!!!!!!!!!!!!!!!!
!
! Dissociate a vector of modes
!
!!!!!!!!!!!!!!!!!

subroutine diss_vec_point_1D(point,N)
use m_derived_types, only : vec_point
implicit none
integer, intent(in) :: N
type(vec_point), intent(out) :: point(N)
! internal
integer :: i

do i=1,N
  nullify(point(i)%w)
enddo

end subroutine diss_vec_point_1D


!!!!!!!!!!!!!!!!!
!
! some simple functions
!
!!!!!!!!!!!!!!!!!

integer function numer_order_param(name)
use m_lattice, only : my_order_parameters
implicit none
character(len=*), intent(in) :: name
! internal
integer :: i

numer_order_param=0
do i=1,size(my_order_parameters)
   if (name.eq.trim(my_order_parameters(i)%name)) numer_order_param=numer_order_param+1
enddo

end function

function find_position_order(N_order_found,name)
use m_lattice, only : my_order_parameters
implicit none
integer, intent(in) :: N_order_found
character(len=*), intent(in) :: name
integer :: find_position_order(3,N_order_found)
! internal
integer :: i,N_orders,N

N_orders=size(my_order_parameters)
N=0
find_position_order=0

do i=1,N_orders
   if (name.eq.trim(my_order_parameters(i)%name)) then
     N=N+1
     find_position_order(1,N_order_found)=i
     find_position_order(2,N_order_found)=my_order_parameters(i)%start
     find_position_order(3,N_order_found)=my_order_parameters(i)%end
   endif
enddo

if (N.ne.N_order_found) stop 'ERROR in find_position_order in module m_operator_pointer_utils'

end function

end module m_operator_pointer_utils
