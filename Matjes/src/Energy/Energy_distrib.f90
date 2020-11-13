module m_energyfield
#if 0
use m_basic_types, only : vec_point, site_ham
use m_derived_types, only : operator_real_order_N,lattice,Op_real_order_N,point_shell_Operator
use m_modes_variables, only : point_shell_mode
use m_Hamiltonian_variables, only : H_vois,shell_Ham_order_N
use m_operator_pointer_utils

! new type for the density distribution
type operator_distrib
   type(operator_real_order_N),allocatable, dimension(:) :: operators
end type operator_distrib

! shell resolved energy tensor
type(operator_distrib), public, protected :: energy_distrib

interface associate_shell_ham
  module procedure associate_shell_ham_1D
end interface associate_shell_ham

private
public :: init_Energy_distrib,get_Energy_distrib
contains

!!!!!!!!!!!!!!!!!!!!!!!
!
! initialize the energy distribution
!
!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_Energy_distrib(my_lattice,tableNN,indexNN)
use m_energy_commons, only : total_hamiltonian
implicit none
type(lattice), intent(in) :: my_lattice
integer, intent(in) :: tableNN(:,:,:,:,:,:) !!tableNN(4,N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),count(my_motif%i_mom)
integer, intent(in) :: indexNN(:)
! input
integer :: size_ham,Nspin,all_size(4)
! slope variables
integer :: i,j
integer,allocatable :: int_ind(:)

all_size=shape(my_lattice%ordpar%l_modes)
Nspin=product(all_size)

size_ham=size(total_hamiltonian%num)
allocate(int_ind(size_ham))
do i=1,size_ham
   int_ind(i)=total_hamiltonian%num(i)
enddo

allocate(energy_distrib%operators(size_ham))

do j=1,size_ham
   call associate_shell_ham(j,energy_distrib%operators(j),total_hamiltonian%shell_num(j),my_lattice,tableNN,int_ind)
enddo

end subroutine init_Energy_distrib







!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine that associates the Hamiltonians for each shells
!
!!!!!!!!!!!!!!!!!!!!!!!

subroutine associate_shell_ham_1D(n_shell,ham_point,ham_target,my_lattice,tableNN,indexNN)
implicit none
type(operator_real_order_N),intent(inout) :: ham_point
type(shell_Ham_order_N),target,intent(in) :: ham_target
type(lattice), intent(in) :: my_lattice
integer, intent(in) :: tableNN(:,:,:,:,:,:) !!tableNN(4,N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),count(my_motif%i_mom)
integer, intent(in) :: indexNN(:),n_shell
! internal variable
integer :: avant,n_atom_shell,Nspin,all_size(4)

all_size=shape(my_lattice%ordpar%l_modes)
Nspin=product(all_size)
n_atom_shell=indexNN(n_shell)

! if the number of atom in the shell is one it means that it is the onsite term.
if (n_atom_shell.eq.1) avant=0

! n_shell is the shell number
if (n_shell.eq.2) then
  avant=0
else
avant=sum(indexNN(2:n_shell-1))
endif

allocate(ham_point%value(n_atom_shell,Nspin))
allocate(ham_point%line(n_atom_shell,Nspin))
ham_point%nline=Nspin
ham_point%ncolumn=Nspin
ham_point%line=0

call associate_pointer(ham_point,ham_target,my_lattice,tableNN,n_atom_shell,avant)

end subroutine associate_shell_ham_1D













subroutine local_energy_pointer_EDestrib(E_int,iomp,i_shell,spin,dim_mode)
use m_dipole_energy
use m_matrix, only : reduce
implicit none
! input
type(vec_point), intent(in) :: spin(:)
integer, intent(in) :: iomp,i_shell,dim_mode
! ouput
real(kind=8), intent(out) :: E_int
! internal
integer :: i,N,j,k
real(kind=8) :: S_int(dim_mode)

N=size(energy_distrib%operators(i_shell)%line(:,iomp))
E_int=0.0d0

do i=1,N

   j=energy_distrib%operators(i_shell)%line(i,iomp)

   call reduce(energy_distrib%operators(i_shell)%value(i,iomp),size(energy_distrib%operators(i_shell)%value(i,iomp)%order_op),S_int,spin(iomp)%w,spin(j)%w,dim_mode)

   E_int=E_int+dot_product( spin(iomp)%w , S_int )

enddo

end subroutine local_energy_pointer_EDestrib






!!!!!!!!!!!!!!!!!!!!!!!
!
! Print the energy density distribution interaction and shell resolve into a file
!
!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_Energy_distrib(tag,spin)
use m_io_files_utils
use m_convert
use m_dipole_energy
use m_dipolar_field, only : i_dip
implicit none
! input
type(vec_point),intent(in) :: spin(:)
integer, intent(in) :: tag
! internals
integer :: iomp,N,io,i_shell
real(kind=8) :: E_int,E_dip
real(kind=8), allocatable :: E_shell(:)
integer :: number_shell,dim_mode
!   name of files
character(len=30) :: fname,rw_format

number_shell=size(energy_distrib%operators)
allocate(E_shell(number_shell))

write(rw_format,'( "(", I4, "(2x,E20.12E3))" )') number_shell+1

E_shell=0.0d0
E_dip=0.0d0
N=size(spin)
dim_mode=size(spin(1)%w)

fname=convert('EnDistrib_',tag,'.dat')
io=open_file_write(fname)

do iomp=1,N

   E_shell=0.0d0

   do i_shell=1,number_shell

      call local_energy_pointer_EDestrib(E_int,iomp,i_shell,spin,dim_mode)
      E_shell(i_shell)=E_shell(i_shell)+E_int

   enddo

   if (i_dip) E_dip=get_dipole_E(iomp)

   write(io,rw_format) (E_shell(i_shell),i_shell=1,number_shell),E_dip
enddo

call close_file(fname,io)

end subroutine get_Energy_distrib
#endif
end module
