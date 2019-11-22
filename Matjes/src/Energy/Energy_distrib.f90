module m_energyfield
use m_derived_types
use m_operator_pointer_utils

! new type for the density distribution
type operator_distrib
   type(operator_real),allocatable, dimension(:) :: operators
   type(point_shell_Operator), allocatable, dimension(:,:) :: E_line
   type(point_shell_mode), allocatable, dimension(:,:) :: mode_E_column
end type operator_distrib

! shell resolved energy tensor
type(operator_distrib), public, protected :: energy_distrib

interface associate_shell_ham
  module procedure associate_shell_ham_1D
end interface associate_shell_ham

private
public :: init_Energy_distrib,get_Energy_distrib,get_Energy_distrib_line
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
integer :: i,j,k
integer,allocatable :: int_ind(:)

all_size=shape(my_lattice%l_modes)
Nspin=product(all_size)
size_ham=size(total_hamiltonian(1)%atom(1)%H,1)


size_ham=size(total_hamiltonian)
allocate(int_ind(size_ham))
do i=1,size_ham
   int_ind(i)=size(total_hamiltonian(i)%atom)
enddo

allocate(energy_distrib%operators(size_ham))

do j=1,size_ham
   call associate_shell_ham(j,energy_distrib%operators(j),total_hamiltonian(j)%atom,my_lattice,tableNN,int_ind)
enddo

!!!!! taking care of the line of the spin matric and column of the energy
allocate(energy_distrib%E_line(Nspin,size_ham),energy_distrib%mode_E_column(Nspin,size_ham))

do i=1,size_ham
   do j=1,Nspin
      allocate(energy_distrib%E_line(j,i)%shell(int_ind(i)),energy_distrib%mode_E_column(j,i)%shell(int_ind(i)))
   enddo
   call dissociate(energy_distrib%E_line(:,i),int_ind(i),Nspin)
   call dissociate(energy_distrib%mode_E_column(:,i),int_ind(i),Nspin)
enddo

end subroutine init_Energy_distrib

!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine that associates the Hamiltonians for each shells
!
!!!!!!!!!!!!!!!!!!!!!!!

subroutine associate_shell_ham_1D(n_shell,ham_point,ham_target,my_lattice,tableNN,indexNN)
implicit none
type(operator_real),intent(inout) :: ham_point
type(site_ham),target,intent(in) :: ham_target(:)
type(lattice), intent(in) :: my_lattice
integer, intent(in) :: tableNN(:,:,:,:,:,:) !!tableNN(4,N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),count(my_motif%i_mom)
integer, intent(in) :: indexNN(:),n_shell
! internal variable
integer :: i,avant,n_atom_shell,Nspin,all_size(4)
integer :: N,M

all_size=shape(my_lattice%l_modes)
Nspin=product(all_size)
M=size(ham_target)
n_atom_shell=indexNN(n_shell)

! if the nu,ber of atom in the shell is one it means that it is the onsite term.
if (n_atom_shell.eq.1) avant=0

! n_shell is the shell number
if (n_shell.eq.2) then
  avant=0
else
  avant=sum(indexNN(2:n_shell-1))
endif

if (M.ne.n_atom_shell) stop 'The Hamiltonian do not have the same size - associate_shell_ham'

allocate(ham_point%value(n_atom_shell,Nspin))
allocate(ham_point%line(n_atom_shell,Nspin))
ham_point%nline=Nspin
ham_point%ncolumn=Nspin
ham_point%line=0

call dissociate(ham_point,Nspin,n_atom_shell)

call associate_pointer(ham_point,ham_target,my_lattice,tableNN,n_atom_shell,avant)

end subroutine associate_shell_ham_1D


!!!!!!!!!!!!!!!!!!!!!!!
!
! associate the line of the Hamiltonians and column of the modes for the energy distribution calculation
!
!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_Energy_distrib_line(spin)
implicit none
! input
type(vec_point),intent(in) :: spin(:)
! internals
integer :: size_line(2),i

!!!!!!!!!!!!!!!
!
! Associate the exchange energy distribution
!
!!!!!!!!!!!!!!!

size_line=shape(energy_distrib%E_line)

do i=1,size_line(2)
   call associate_pointer(energy_distrib%mode_E_column(:,i),spin,energy_distrib%E_line(:,i),energy_distrib%operators(i))
enddo

end subroutine get_Energy_distrib_line





!!!!!!!!!!!!!!!!!!!!!!!
!
! Print the energy density distribution interaction and shell resolve into a file
!
!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_Energy_distrib(tag,spin,h_int)
use m_io_files_utils
use m_convert
!use  m_total_energy, only : total_energy
use m_local_energy, only : local_energy_pointer
implicit none
! input
type(vec_point),intent(in) :: spin(:)
integer, intent(in) :: tag
real(kind=8), intent(in) :: h_int(:)
! internals
integer :: i,N,io,j
real(kind=8) :: E_int
real(kind=8), allocatable :: E_shell(:)
integer :: number_shell
!   name of files
character(len=30) :: fname,rw_format

number_shell=size(energy_distrib%operators)
allocate(E_shell(number_shell))

write(rw_format,'( "(", I4, "(2x,f20.15))" )') number_shell

E_shell=0.0d0
N=size(spin)

fname=convert('EnDistrib_',tag,'.dat')
io=open_file_write(fname)

do i=1,N

   E_shell=0.0d0
   do j=1,number_shell

      call local_energy_pointer(E_int,i,energy_distrib%mode_E_column(i,j),energy_distrib%E_line(i,j))
      E_shell(j)=E_shell(j)+E_int

   enddo

   write(io,rw_format) E_shell
enddo
#ifdef CPP_BRUTDIP
!          if (i_dip) E_dip=dipole(i_x,i_y,i_z,i_m,spin,shape_spin,masque,my_lattice%boundary)
#else
!          if (i_dip) E_dip=fftdip(i_x,i_y,i_z,i_m)
#endif

!          Write(70,'(8(E20.12E3,2x))') E_DM,E_ani,E_z,E_4,E_biq,E_dip

call close_file(fname,io)

end subroutine get_Energy_distrib

end module
