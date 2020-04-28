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
   module procedure convoluate_Op_sym_file
end interface

private
public :: get_diagonal_Op,get_symmetric_Op,get_Op_in_Op,convoluate_Op_SOC_vector,get_borders,convoluate_Op_sym

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
subroutine convoluate_Op_sym_file(k,H_in,H_out,x_start,x_end,y_start,y_end)
use m_basic_types, only : symop
use m_grp_sym
use m_matrix
implicit none
integer, intent(in) :: k,x_start,x_end,y_start,y_end
real(kind=8), intent(in) :: H_in(:,:)
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
! subroutine that print stupid error message
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine error(string)
implicit none
character(len=*), intent(in) :: string

write(6,'(2a)') 'error in the subroutine  ',string
stop

end subroutine







! dumy functions
subroutine get_borders(name_int,x_start,x_end,name_ext,y_start,y_end,my_order_parameters)
use m_derived_types, only : order_parameter
implicit none
character(len=*), intent(in) :: name_int,name_ext
type(order_parameter), intent(in) :: my_order_parameters(:)
integer, intent(out) :: x_start,x_end,y_start,y_end
! internal
integer :: n_mode_int,n_mode_ext
integer, allocatable :: positions_int(:,:), positions_ext(:,:)

n_mode_int=find_n_mode(name_int,my_order_parameters)
allocate(positions_int(2,n_mode_int))
positions_int=0
call find_position(name_int,my_order_parameters,positions_int)

x_start=minval(positions_int)
x_end=maxval(positions_int)

if (((x_end-x_start+1)/3).ne.size(positions_int,2)) stop 'ERROR in couplage_ME - data not contiguous'

n_mode_ext=find_n_mode(name_ext,my_order_parameters)
allocate(positions_ext(2,n_mode_ext))
positions_ext=0
call find_position(name_ext,my_order_parameters,positions_ext)
! check the data

y_start=minval(positions_ext)
y_end=maxval(positions_ext)

if (((y_end-y_start+1)/3).ne.size(positions_ext,2)) stop 'ERROR in couplage_ME - data not contiguous'

end subroutine


! find the position of the order parameters in the lattice
subroutine find_position(name,my_order_parameters,position)
use m_derived_types, only : order_parameter
implicit none
character(len=*), intent(in) :: name
type(order_parameter), intent(in) :: my_order_parameters(:)
integer, intent(inout) :: position(:,:)
! internal
integer :: n_mode,i,n_mode_int

n_mode=size(position,2)
n_mode_int=0

do i=1,size(my_order_parameters)
   if (name.eq.trim(my_order_parameters(i)%name)) then
      n_mode_int=n_mode_int+1
      position(1,n_mode_int)=my_order_parameters(i)%start
      position(2,n_mode_int)=my_order_parameters(i)%end
   endif
enddo

if (n_mode_int.ne.n_mode) stop 'ERROR in find_position - more mode read than found'

end subroutine

! find the occurence of name in an array of order parameters
function find_n_mode(name,my_order_parameters)
use m_derived_types, only : order_parameter
implicit none
character(len=*), intent(in) :: name
type(order_parameter), intent(in) :: my_order_parameters(:)
integer :: find_n_mode
! internal
integer :: i

find_n_mode=0
! find the number of occurence of name in the array of my_order_parameters
do i=1,size(my_order_parameters)
   if (name.eq.trim(my_order_parameters(i)%name)) find_n_mode=find_n_mode+1
enddo

end function

end module m_symmetry_operators
