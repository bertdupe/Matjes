module m_summer_exp
use m_lattice, only : my_order_parameters
use m_Hamiltonian_variables, only : coeff_ham_inter_spec
use m_derived_types
use m_derived_types, only : lattice
use m_convert

! do we turn on the interaction
logical :: i_temperature_Stra=.False.

! coupling term between T and the Js at d=0nm (onsite term)
real(kind=8) :: T_coeff=0.0d0
real(kind=8) :: X_coeff=0.0d0
integer :: n_coeff=1
real(kind=8), allocatable, dimension(:) :: T_coeff_shell

type(coeff_ham_inter_spec), target, public, protected :: temperature_strasbourg






private
public :: get_coeff_TStra,get_Temperature_H

contains

!
! Initilization routine
!

subroutine get_coeff_TStra(my_lattice,my_motif)
use m_io_utils
use m_io_files_utils
use m_table_dist
use m_tbessy
use m_tbessj
use m_tbessk
implicit none
type(lattice), intent(in) :: my_lattice
type(t_cell), intent(in) :: my_motif
! internal
integer :: io_input,i
real(kind=8), allocatable, dimension(:,:) :: tabledist
character(len=50) :: form
real(kind=8) :: alpha

temperature_strasbourg%i_exist=i_temperature_Stra

io_input=open_file_read('input')

call get_parameter(io_input,'input','T_stra_on',i_temperature_Stra)
call get_parameter(io_input,'input','T_coeff',T_coeff)
call get_parameter(io_input,'input','X_coeff',X_coeff)
call get_parameter(io_input,'input','n_T_coeff',n_coeff)

call close_file('input',io_input)

if (.not.i_temperature_Stra) return

if (abs(X_coeff).lt.1.0d-8) stop 'X_coeff must be different from 0'

allocate(tabledist(n_coeff,1),T_coeff_shell(n_coeff+1))
tabledist=0.0d0
T_coeff_shell=0.0d0
T_coeff_shell(1)=0.0d0

call get_table_of_distance(my_lattice%areal,n_coeff,my_lattice%world,my_motif,tabledist)
do i=1,n_coeff
  alpha=X_coeff*tabledist(i,1)
  T_coeff_shell(i+1)=T_coeff*(BESSJ(0,alpha)*BESSK(1,alpha)-BESSJ(1,alpha)*BESSK(0,alpha))
enddo

form=convert('(',n_coeff+1,'(f14.8,2x))')
write(6,'(a)') ''
write(6,'(a)') 'Temperature correction of the Hamiltonian'
write(6,form) T_coeff_shell
write(6,'(a)') ''

end subroutine

!
! update of the energy
!

subroutine get_Temperature_H(dim_ham)
use m_symmetry_operators
implicit none
integer, intent(in) :: dim_ham
! internal
integer :: n_shell
integer :: i,j
! temperature
integer :: x_start,x_end
! temperature
integer :: y_start,y_end
character(len=50) :: form

if (.not.i_temperature_Stra) return

n_shell=size(T_coeff_shell)
temperature_strasbourg%name='temp-stra'
temperature_strasbourg%c_ham=1.0d0
temperature_strasbourg%N_shell=n_shell
temperature_strasbourg%order=2
temperature_strasbourg%i_exist=i_temperature_Stra

allocate(temperature_strasbourg%ham(n_shell))
do i=1,n_shell
   allocate(temperature_strasbourg%ham(i)%H(dim_ham,dim_ham))
   temperature_strasbourg%ham(i)%H=0.0d0
enddo

call get_borders('magnetic',x_start,x_end,'magnetic',y_start,y_end,my_order_parameters)

do i=1,n_shell
  call get_diagonal_Op(temperature_strasbourg%ham(i)%H,T_coeff_shell(i),temperature_strasbourg%c_ham,x_start,x_end)
enddo


form=convert('(',dim_ham,'(f12.8,2x))')
write(6,'(a)') ''
write(6,'(a)') 'Temperature tensor of order 2'
do j=1,n_shell
  write(6,'(a,I8)') 'shell ', j-1
  do i=1,dim_ham
    write(6,form) temperature_strasbourg%ham(j)%H(:,i)
  enddo
enddo
write(6,'(a)') ''

end subroutine get_Temperature_H


end module m_summer_exp
