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

private
public :: get_diagonal_Op,get_symmetric_Op,get_Op_in_Op

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
