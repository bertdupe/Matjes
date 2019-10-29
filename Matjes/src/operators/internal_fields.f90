module m_internal_fields_commons
use m_derived_types, only : operator_real,shell_Ham,point_shell_Operator,vec_point,point_shell_mode
use m_operator_pointer_utils
! coefficients of the internal magnetic field
type(shell_Ham), public, protected, target, save, allocatable, dimension(:) :: B_int
!type(shell_Ham), public, protected, target, save, allocatable, dimension(:) :: total_hamiltonian

! total effective field tensor
type(operator_real),public,protected,save :: B_total

private
public :: associate_internal_Beff,get_B_line

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
use m_derived_types, only : Coeff_Ham,lattice
use m_energy_commons, only : energy,total_hamiltonian
use m_constants, only : identity
use m_get_position
implicit none
type(lattice), intent(in) :: my_lattice
integer, intent(in) :: tableNN(:,:,:,:,:,:) !!tableNN(4,N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),count(my_motif%i_mom)
integer, intent(in) :: indexNN(:)
! internal
!type(operator_real),save :: energy
! slope of the sums
integer :: i,j,N_coeff_DMI,N_all_vois,N_all_shell,N_coeff_ani,N_coeff_Exch,size_energy(2)
integer :: Nl,Nc,dim_ham

Nl=energy%nline
Nc=energy%ncolumn

!Hamiltonian%total_shell(i)%atom(l)%H(2,1)
! allocate the variables
N_all_shell=size(total_hamiltonian)
dim_ham=my_lattice%dim_mode

allocate(B_int(N_all_shell))
do i=1,N_all_shell
   N_all_vois=size(total_hamiltonian(i)%atom)
   allocate(B_int(i)%atom(N_all_vois))
   do j=1,N_all_vois
      allocate(B_int(i)%atom(j)%H(dim_ham,dim_ham))
      B_int(i)%atom(j)%H=0.0d0
   enddo
enddo

do i=1,N_all_shell

   N_all_vois=size(total_hamiltonian(i)%atom)
   do j=1,N_all_vois

      B_int(i)%atom(j)%H=-2.0d0*total_hamiltonian(i)%atom(j)%H

   enddo

enddo

size_energy=shape(energy%value)
allocate(B_total%value(size_energy(1),size_energy(2)))
allocate(B_total%line(size_energy(1),size_energy(2)))
B_total%nline=Nl
B_total%ncolumn=Nc
B_total%line=0

! first nullify all pointer
call dissociate(B_total%value,size_energy(1),size_energy(2))

call associate_pointer(B_total,B_int,my_lattice,tableNN,indexNN)

#ifdef CPP_DEBUG
do i=1,N_all_shell
   write(*,*) 'shell', i

   do j=1,N_all_vois
      if (all( abs(Hamiltonian%total_shell(i)%atom(j)%H(:,:)).lt.1.0d-8 )) cycle

      write(6,*) B_int%total_shell(i)%atom(j)%H(1,:)
      write(6,*) B_int%total_shell(i)%atom(j)%H(2,:)
      write(6,*) B_int%total_shell(i)%atom(j)%H(3,:)
      pause
   enddo
enddo
#endif

end subroutine associate_internal_Beff

subroutine get_B_line(B_line,mode_B_column,spin)
use m_derived_types, only : point_shell_Operator,point_shell_mode,vec_point
use m_operator_pointer_utils
implicit none
type(point_shell_Operator), intent(inout) :: B_line(:)
type(point_shell_mode),intent(inout) :: mode_B_column(:)
type(vec_point),target,intent(in) :: spin(:)
! internal variables
integer :: shape_B_total(2),i

shape_B_total=shape(B_total%value)

do i=1,shape_B_total(2)
   allocate(B_line(i)%shell(shape_B_total(1)))
   allocate(mode_B_column(i)%shell(shape_B_total(1)))
enddo
call dissociate(B_line,shape_B_total(1),shape_B_total(2))
call dissociate(mode_B_column,shape_B_total(1),shape_B_total(2))


call associate_pointer(mode_B_column,spin,B_line,B_total)

end subroutine get_B_line

end module m_internal_fields_commons
