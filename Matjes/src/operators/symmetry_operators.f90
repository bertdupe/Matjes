module m_symmetry_operators

interface get_diagonal_Op
   module procedure get_diagonal_Op_real_3D_c_2D,get_diagonal_Op_real_3D_c_0D
end interface

interface get_symmetric_Op
   module procedure get_bloc_diag_symmetric_Op_real_2D
end interface

interface get_Op_in_Op
   module procedure get_Op_in_Op_Xshell
end interface

interface convoluate_Op_SOC_vector
   module procedure convoluate_Op_2D_SOC_vector_1D
end interface

private
public :: get_diagonal_Op,get_symmetric_Op,get_Op_in_Op,convoluate_Op_SOC_vector

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
subroutine get_diagonal_Op_real_3D_c_0D(Op,coeff,constant,start,end)
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
subroutine get_Op_in_Op_Xshell(Op_big,Op_small,upper_left_corner_x,upper_left_corner_y,start,end)
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
subroutine convoluate_Op_2D_SOC_vector_1D(D,Op_DMI,Op_DMI_Nshell_local)
implicit none
real(kind=8), intent(in) :: D(:,:),Op_DMI(:,:,:)
real(kind=8), intent(out) :: Op_DMI_Nshell_local(:,:,:,:)
! internal variable
integer :: shape_Op_dmi(3),dim_D(2),i,j

shape_Op_dmi=shape(Op_DMI)
dim_D=shape(D)

Op_DMI_Nshell_local=0.0d0

do i=1,shape_Op_dmi(3)  ! loop over the shells
   do j=1,dim_D(1)    ! loop of each atom in one shell
      Op_DMI_Nshell_local(1,2,j,i)=Op_DMI(1,2,i)*D(j,3)
      Op_DMI_Nshell_local(1,3,j,i)=Op_DMI(1,3,i)*D(j,2)
      Op_DMI_Nshell_local(2,3,j,i)=Op_DMI(2,3,i)*D(j,1)

      Op_DMI_Nshell_local(2,1,j,i)=Op_DMI(2,1,i)*D(j,3)
      Op_DMI_Nshell_local(3,1,j,i)=Op_DMI(3,1,i)*D(j,2)
      Op_DMI_Nshell_local(3,2,j,i)=Op_DMI(3,2,i)*D(j,1)
   enddo
enddo

end subroutine convoluate_Op_2D_SOC_vector_1D





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! get the off diagonal symmetric part of the Hamiltonian where you have 1 shell involved
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_bloc_diag_symmetric_Op_real_2D(Op,constant,upper_left_corner_x,upper_left_corner_y,start,end)
implicit none
real(kind=8), intent(inout) :: Op(:,:)
real(kind=8), intent(in) :: constant
integer, intent(in) :: start,end,upper_left_corner_x,upper_left_corner_y
! internal
integer :: shape_Op(2)
integer :: j,index_mode_x,index_mode_y

shape_Op=shape(Op)

index_mode_x=upper_left_corner_x-1
index_mode_y=upper_left_corner_y

do j=1,end-start+1
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
