module m_internal_fields_commons
use m_basic_types, only : vec_point
use m_derived_types, only : operator_real_order_N,point_shell_Operator
use m_Hamiltonian_variables, only : shell_Ham,H_vois
use m_modes_variables, only : point_shell_mode
use m_operator_pointer_utils
! coefficients of the internal magnetic field
type(H_vois), public, protected, target, save :: B_int

! total effective field tensor
type(operator_real_order_N),public,protected,save :: B_total

private
public :: associate_internal_Beff

contains

!!!!!!!!!!!!!!!!!!!!!!!
!  Routine that prepares the LOCAL effective field based on the Hamiltonian in
!  energy_commons.f90
! The form of the effective field is the same as the Hamiltonian
!
!  ( B_00  ...           )  (Site_0)
!  ( B_10  B_11          )  (Site_1)
!  (  .          .       )  (Site_2)
!  (  .             .    )  (...   )
!  (  .               .  )  (...   )
!
!!!!!!!!!!!!!!!!!!!!!!!

subroutine associate_internal_Beff(my_lattice,tableNN,indexNN)
use m_derived_types, only : lattice
use m_Hamiltonian_variables, only : Coeff_Ham
use m_energy_commons, only : energy,total_hamiltonian
use m_constants, only : identity
use m_get_position
use m_convert
implicit none
type(lattice), intent(in) :: my_lattice
integer, intent(in) :: tableNN(:,:,:,:,:,:) !!tableNN(4,N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),count(my_motif%i_mom)
integer, intent(in) :: indexNN(:)
! internal
integer :: i,j,N_coeff_DMI,N_all_vois,N_all_shell,N_coeff_ani,N_coeff_Exch,size_energy(2),k,l
integer :: Nl,Nc,dim_ham,N_order,i_order,order
character(len=50) :: form

Nl=energy%nline
Nc=energy%ncolumn

! allocate the variables
N_all_shell=size(total_hamiltonian%shell_num)
dim_ham=my_lattice%dim_mode

allocate(B_int%shell_num(N_all_shell),B_int%num(N_all_shell))
B_int%num=total_hamiltonian%num
do i=1,N_all_shell
   N_all_vois=B_int%num(i)
   N_order=size(total_hamiltonian%shell_num(i)%num)
   allocate(B_int%shell_num(i)%order(N_order),B_int%shell_num(i)%num(N_order))
   B_int%shell_num(i)%num=total_hamiltonian%shell_num(i)%num

   do i_order=1,N_order

     allocate(B_int%shell_num(i)%order(i_order)%atom(N_all_vois))
     order=B_int%shell_num(i)%num(i_order)
     B_int%shell_num(i)%order(i_order)%line=dim_ham**( order-1 )
     B_int%shell_num(i)%order(i_order)%column=dim_ham

     do j=1,N_all_vois

       allocate(B_int%shell_num(i)%order(i_order)%atom(j)%H( dim_ham , dim_ham**( order-1 ) ))
       B_int%shell_num(i)%order(i_order)%atom(j)%H=0.0d0

     enddo
   enddo
enddo

do i=1,N_all_shell

   N_all_vois=B_int%num(i)
   N_order=size(B_int%shell_num(i)%num)
   do i_order=1,N_order
     do j=1,N_all_vois

       B_int%shell_num(i)%order(i_order)%atom(j)%H=-2.0d0*total_hamiltonian%shell_num(i)%order(i_order)%atom(j)%H

     enddo
   enddo
enddo

size_energy=shape(energy%value)
allocate(B_total%value(size_energy(1),size_energy(2)))
allocate(B_total%line(size_energy(1),size_energy(2)))
B_total%nline=Nl
B_total%ncolumn=Nc
B_total%line=0

call associate_pointer(B_total,B_int,my_lattice,tableNN,indexNN)

form=convert('(',dim_ham,'(2x,f12.8))')
write(6,'(a)') ''
write(6,'(a)') 'Total crystal field B_eff'
do i=1,size(B_int%num)
  write(6,'(a,2x,I3)') 'B_eff for shell ',i

  do j=1,size(B_int%shell_num(i)%num)
    write(6,'(a,2x,I3)') 'B_eff of order ',j

    do k=1,B_int%num(i)

        write(6,'(a,2x,I3)') 'Bond', k

      do l=1,B_int%shell_num(i)%order(j)%line
        write(6,form) B_int%shell_num(i)%order(j)%atom(k)%H(:,l)
      enddo

    enddo
    write(6,'(a)') ''

  enddo
  write(6,'(a)') ''
enddo

end subroutine associate_internal_Beff

!subroutine get_B_line(B_line,mode_B_column,spin)
!use m_operator_pointer_utils
!implicit none
!type(point_shell_Operator), intent(inout) :: B_line(:)
!type(point_shell_mode),intent(inout) :: mode_B_column(:)
!type(vec_point),target,intent(in) :: spin(:)
!! internal variables
!integer :: shape_B_total(2),i

!shape_B_total=shape(B_total%value)
!
!do i=1,shape_B_total(2)
!   allocate(B_line(i)%shell(shape_B_total(1)))
!   allocate(mode_B_column(i)%shell(shape_B_total(1)))
!enddo
!call dissociate(B_line,shape_B_total(1),shape_B_total(2))
!call dissociate(mode_B_column,shape_B_total(1),shape_B_total(2))
!
!
!call associate_pointer(mode_B_column,spin,B_line,B_total)

!end subroutine get_B_line

end module m_internal_fields_commons
