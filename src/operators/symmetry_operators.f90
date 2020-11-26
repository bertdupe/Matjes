module m_symmetry_operators

interface get_diagonal_Op
   module procedure get_diagonal_Op_real_3D_c_2D,get_diagonal_Op_real_3D_c_1D,get_diagonal_Op_real_2D_c_0D
end interface

interface get_symmetric_Op
   module procedure get_bloc_diag_symmetric_Op_real_2D
end interface

interface get_Op_in_Op
   module procedure get_Op_3D_in_Op_3D_Xshell,get_Op_2D_in_Op_2D
end interface

interface convoluate_Op_SOC_vector
   module procedure convoluate_Op_2D_SOC_vector_1D
end interface

interface convoluate_Op_sym
   module procedure convoluate_Op_sym_file,convoluate_Op_sym_bond
end interface

private
public :: get_diagonal_Op,get_symmetric_Op,get_Op_in_Op,convoluate_Op_SOC_vector,convoluate_Op_sym

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! get the diagonal part of the Hamiltonian where you have more than 1 shell involved and the coefficients are fields
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_diagonal_Op_real_3D_c_2D(Op,coeff,constant,start,end)
implicit none
real(kind=8), intent(inout) :: Op(:,:,:)
real(kind=8), intent(in) :: coeff(:,:)
real(kind=8), intent(in) :: constant
integer, intent(in) :: start,end
! internal
integer :: shape_Op(3), shape_coeff(2)
integer :: i,j,index_mode

shape_Op=shape(Op)
shape_coeff=shape(coeff)

if (shape_coeff(1).ne.(end-start+1)) call error('get_diagonal_Op_real_3D - the mode dimension does not match the Operator')
if (shape_coeff(2).ne.shape_Op(3)) call error('get_diagonal_Op_real_3D - number of shell incorrect')

do i=1,shape_Op(3)
   index_mode=start
   do j=1,end-start+1
      Op(index_mode,index_mode,i)=coeff(j,i)*constant
      index_mode=index_mode+1
   enddo
enddo

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! get the diagonal part of the Hamiltonian where you have more than 1 shell involved and the coefficients are constants
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_diagonal_Op_real_2D_c_0D(Op,coeff,constant,start,end)
implicit none
real(kind=8), intent(inout) :: Op(:,:)
real(kind=8), intent(in) :: coeff
real(kind=8), intent(in) :: constant
integer, intent(in) :: start,end
! internal
integer :: shape_Op(2)
integer :: j,index_mode

shape_Op=shape(Op)
index_mode=start

do j=1,end-start+1
   Op(index_mode,index_mode)=coeff*constant
   index_mode=index_mode+1
enddo

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! get the diagonal part of the Hamiltonian where you have more than 1 shell involved and the coefficients are constants
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_diagonal_Op_real_3D_c_1D(Op,coeff,constant,start,end)
implicit none
real(kind=8), intent(inout) :: Op(:,:,:)
real(kind=8), intent(in) :: coeff(:)
real(kind=8), intent(in) :: constant
integer, intent(in) :: start,end
! internal
integer :: shape_Op(3), shape_coeff(1)
integer :: i,j,index_mode

shape_Op=shape(Op)
shape_coeff=shape(coeff)

if (shape_coeff(1).ne.shape_Op(3)) call error('get_diagonal_Op_real_3D - number of shell incorrect')

do i=1,shape_Op(3)
   index_mode=start
   do j=1,end-start+1
      Op(index_mode,index_mode,i)=coeff(i)*constant
      index_mode=index_mode+1
   enddo
enddo

end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! get a small Hamiltonian into a bigger one
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_Op_2D_in_Op_2D(Op_big,Op_small,upper_left_corner_x,down_left_corner_x,upper_left_corner_y,down_left_corner_y)
implicit none
real(kind=8), intent(inout) :: Op_big(:,:)
real(kind=8), intent(in) :: Op_small(:,:)
integer, intent(in) :: upper_left_corner_y,down_left_corner_y,upper_left_corner_x,down_left_corner_x
! internal
integer :: shape_Op_big(2), shape_Op_small(2)
integer :: j,k,index_mode_x,index_mode_y

shape_Op_big=shape(Op_big)
shape_Op_small=shape(Op_small)

if ((shape_Op_small(1).ne.(down_left_corner_y-upper_left_corner_y+1)).or.(shape_Op_small(2).ne.(down_left_corner_y-upper_left_corner_y+1))) call error('get_Op_in_Op - size Op incorrect')

index_mode_y=upper_left_corner_y-1

do j=1,shape_Op_small(2)
   index_mode_y=index_mode_y+1
   index_mode_x=upper_left_corner_x-1

   do k=1,shape_Op_small(1)
      index_mode_x=index_mode_x+1
      Op_big(index_mode_x,index_mode_y)=Op_small(k,j)+Op_big(index_mode_x,index_mode_y)
   enddo

enddo

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! get a small Hamiltonian into a bigger one
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_Op_3D_in_Op_3D_Xshell(Op_big,Op_small,upper_left_corner_x,upper_left_corner_y,start,end)
implicit none
real(kind=8), intent(inout) :: Op_big(:,:,:)
real(kind=8), intent(in) :: Op_small(:,:,:)
integer, intent(in) :: start,end,upper_left_corner_x,upper_left_corner_y
! internal
integer :: shape_Op_big(3), shape_Op_small(3)
integer :: i,j,k,index_mode_x,index_mode_y

shape_Op_big=shape(Op_big)
shape_Op_small=shape(Op_small)

if (shape_Op_big(3).ne.shape_Op_small(3)) call error('get_Op_in_Op_Xshell - number of shell incorrect')
if ((shape_Op_small(1).ne.(end-start+1)).or.(shape_Op_small(2).ne.(end-start+1))) call error('get_Op_in_Op_Xshell - size Op incorrect')

do i=1,shape_Op_big(3)

   index_mode_y=upper_left_corner_y-1

   do j=1,shape_Op_small(2)
     index_mode_y=index_mode_y+1
     index_mode_x=upper_left_corner_x-1

      do k=1,shape_Op_small(1)
        index_mode_x=index_mode_x+1
        Op_big(index_mode_x,index_mode_y,i)=Op_small(k,j,i)
      enddo

   enddo
enddo

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! convolutate the Hamiltonian with the D vector for the DMI rules
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine convoluate_Op_2D_SOC_vector_1D(D,Op_DMI,H)
use m_derived_types, only : site_Ham
implicit none
real(kind=8), intent(in) :: D(:),Op_DMI(:,:)
real(kind=8), intent(inout) :: H(:,:)
! internal variable

H(1,2)=Op_DMI(1,2)*D(3)
H(1,3)=Op_DMI(1,3)*D(2)
H(2,3)=Op_DMI(2,3)*D(1)

H(2,1)=Op_DMI(2,1)*D(3)
H(3,1)=Op_DMI(3,1)*D(2)
H(3,2)=Op_DMI(3,2)*D(1)

H(1,1)=Op_DMI(1,1)
H(2,2)=Op_DMI(2,2)
H(3,3)=Op_DMI(3,3)

end subroutine convoluate_Op_2D_SOC_vector_1D





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! get the off diagonal symmetric part of the Hamiltonian where you have 1 shell involved
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_bloc_diag_symmetric_Op_real_2D(Op,constant,upper_corner_x,upper_corner_y,down_corner_x,down_corner_y)
implicit none
real(kind=8), intent(inout) :: Op(:,:)
real(kind=8), intent(in) :: constant
integer, intent(in) :: down_corner_x,down_corner_y,upper_corner_x,upper_corner_y
! internal
integer :: shape_Op(2)
integer :: j,index_mode_x,index_mode_y

shape_Op=shape(Op)

index_mode_x=upper_corner_x-1
index_mode_y=upper_corner_y-1

if ((down_corner_x-upper_corner_x+1).ne.(down_corner_y-upper_corner_y+1)) stop 'Error in symmetry_operation - get_bloc_diag_symmetric_Op_real_2D'

do j=1,down_corner_x-upper_corner_x+1
   index_mode_x=index_mode_x+1
   index_mode_y=index_mode_y+1
   Op(index_mode_x,index_mode_y)=constant
   Op(index_mode_y,index_mode_x)=constant
enddo

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! convolute a matrix with the k-symmetry operator from file
! This will not work with 2 atoms in one unit cell
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine convoluate_Op_sym_file(k,H_in,H_out,x_start,x_end,y_start,y_end,file)
use m_basic_types, only : symop
use m_grp_sym
use m_invert
implicit none
integer, intent(in) :: k,x_start,x_end,y_start,y_end
real(kind=8), intent(in) :: H_in(:,:)
character(len=*), intent(in) :: file
real(kind=8), intent(inout) :: H_out(:,:)
! internal
integer :: i,j,shape_H(2),l,m,n,o
real(kind=8) :: H_in_local(3,3),H_dum(3,3)
type(symop), allocatable :: symmetries(:), invert_symmetries(:)
real(kind=8), allocatable :: test(:,:,:),test_out(:,:,:),sym_mat(:,:)
integer :: n_sym

H_out=0.0d0

n_sym=get_num_sym_file()
allocate(symmetries(n_sym),invert_symmetries(n_sym))
call read_symmetries(n_sym,symmetries)

do i=1,n_sym
  invert_symmetries(i)%name=symmetries(i)%name
  call invert(symmetries(i)%mat,invert_symmetries(i)%mat,3)
enddo

shape_H=shape(H_in)
allocate(test(shape_H(1),shape_H(1),shape_H(1)),test_out(shape_H(1),shape_H(1),shape_H(1)),sym_mat(shape_H(1),shape_H(1)))
test=reshape(H_in,(/shape_H(1),shape_H(1),shape_H(1)/))

sym_mat=0.0d0
!
! first put ones on the diagonals
!
do i=1,shape_H(1)/3
  sym_mat( 3*(i-1)+1 : 3*i ,3*(i-1)+1 : 3*i )=symmetries(k)%mat
enddo

test_out=0.0d0
do i=1,shape_H(1)
  do j=1,shape_H(1)
    do o=1,shape_H(1)

    do l=1,shape_H(1)
      do m=1,shape_H(1)
        do n=1,shape_H(1)

     test_out(i,j,o)=test_out(i,j,o)+sym_mat(i,l)*sym_mat(j,m)*sym_mat(o,n)*test(l,m,n)

        enddo
      enddo
    enddo

    enddo
  enddo
enddo

H_out=reshape(test_out,(/shape_H(1),shape_H(2)/))

end subroutine convoluate_Op_sym_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! convolute a matrix with the rotation operator for the k-bond
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine convoluate_Op_sym_bond(i_shell,i_bond,H_in,H_out,x_start,x_end,y_start,y_end)
use m_basic_types, only : symop
use m_arrange_neigh, only : get_rotation
implicit none
integer, intent(in) :: i_bond,x_start,x_end,y_start,y_end,i_shell
real(kind=8), intent(in) :: H_in(:,:)
real(kind=8), intent(inout) :: H_out(:,:)
! internal
integer :: i,j,shape_H(2),l,m,n,o
real(kind=8) :: H_dum(3,3)
real(kind=8), allocatable :: test(:,:,:),test_out(:,:,:),sym_mat(:,:)
integer :: n_sym

H_out=0.0d0
H_dum=0.0d0
call get_rotation(i_shell,i_bond,H_dum)

shape_H=shape(H_in)
allocate(test(shape_H(1),shape_H(1),shape_H(1)),test_out(shape_H(1),shape_H(1),shape_H(1)),sym_mat(shape_H(1),shape_H(1)))
test=reshape(H_in,(/shape_H(1),shape_H(1),shape_H(1)/))

sym_mat=0.0d0
!
! first put ones on the diagonals
!
do i=1,shape_H(1)/3
  sym_mat( 3*(i-1)+1 : 3*i ,3*(i-1)+1 : 3*i )=H_dum
enddo

test_out=0.0d0
do i=1,shape_H(1)
  do j=1,shape_H(1)
    do o=1,shape_H(1)

    do l=1,shape_H(1)
      do m=1,shape_H(1)
        do n=1,shape_H(1)

     test_out(i,j,o)=test_out(i,j,o)+sym_mat(i,l)*sym_mat(j,m)*sym_mat(o,n)*test(l,m,n)

        enddo
      enddo
    enddo

    enddo
  enddo
enddo

H_out=reshape(test_out,(/shape_H(1),shape_H(2)/))

end subroutine convoluate_Op_sym_bond














!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine that print stupid error message
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine error(string)
implicit none
character(len=*), intent(in) :: string

write(6,'(2a)') 'error in the subroutine  ',string
stop

end subroutine

end module m_symmetry_operators
