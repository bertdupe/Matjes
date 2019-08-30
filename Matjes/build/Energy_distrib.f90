module m_energyfield
use m_derived_types
use m_operator_pointer_utils

! new type for the density distribution
type operator_distrib
   type(operator_real),allocatable, dimension(:) :: operators
   type(point_shell_Operator), allocatable, dimension(:,:) :: E_line
   type(point_shell_mode), allocatable, dimension(:,:) :: mode_E_column
end type operator_distrib

! Exchange energy tensor
type(operator_distrib) :: exchange
! DMI energy tensor
type(operator_distrib) :: DMI
! anisotropy energy tensor
type(operator_distrib) :: anisotropy
! Zeeman energy tensor
type(operator_real) :: Zeeman

!type(point_shell_Operator), allocatable, dimension(:) :: E_line
!type(point_shell_mode), allocatable, dimension(:) :: mode_E_column

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
use m_energy_commons, only : Hamiltonian
implicit none
type(lattice), intent(in) :: my_lattice
integer, intent(in) :: tableNN(:,:,:,:,:,:) !!tableNN(4,N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),count(my_motif%i_mom)
integer, intent(in) :: indexNN(:)
! input
integer :: size_ham(3),Nspin,all_size(4)
! slope variables
integer :: i,j

all_size=shape(my_lattice%l_modes)
Nspin=product(all_size)

!!!!!!!!!!!!!!!
!
! Exchange energy distribution
!
!!!!!!!!!!!!!!!
size_ham=shape(Hamiltonian%exchange)
allocate(exchange%operators(size_ham(3)))

call associate_shell_ham(exchange%operators,size_ham(3),Hamiltonian%exchange,my_lattice,tableNN,indexNN)

!!!!! taking care of the line of the spin matric and column of the energy
allocate(exchange%E_line(Nspin,size_ham(3)),exchange%mode_E_column(Nspin,size_ham(3)))

do i=1,size_ham(3)
   do j=1,Nspin
      allocate(exchange%E_line(j,i)%shell(indexNN(i)),exchange%mode_E_column(j,i)%shell(indexNN(i)))
   enddo
   call dissociate(exchange%E_line(:,i),indexNN(i),Nspin)
   call dissociate(exchange%mode_E_column(:,i),indexNN(i),Nspin)
enddo

!!!!!!!!!!!!!!!
!
! DMI energy distribution
!
!!!!!!!!!!!!!!!

size_ham=shape(Hamiltonian%DMI)
allocate(DMI%operators(size_ham(3)))

call associate_shell_ham(DMI%operators(:),size_ham(3),Hamiltonian%DMI,my_lattice,tableNN,indexNN)

!!!!! taking care of the line of the spin matric and column of the energy
allocate(DMI%E_line(Nspin,size_ham(3)),DMI%mode_E_column(Nspin,size_ham(3)))

do i=1,size_ham(3)
   do j=1,Nspin
      allocate(DMI%E_line(j,i)%shell(indexNN(i)),DMI%mode_E_column(j,i)%shell(indexNN(i)))
   enddo
   call dissociate(DMI%E_line(:,i),indexNN(i),Nspin)
   call dissociate(DMI%mode_E_column(:,i),indexNN(i),Nspin)
enddo

!!!!!!!!!!!!!!!
!
! anisotropy energy distribution
!
!!!!!!!!!!!!!!!
allocate(anisotropy%operators(size_ham(3)))

call associate_shell_ham(anisotropy%operators(:),size_ham(3),Hamiltonian%ani,my_lattice,tableNN,indexNN)

!!!!! taking care of the line of the spin matric and column of the energy
allocate(anisotropy%E_line(Nspin,size_ham(3)),anisotropy%mode_E_column(Nspin,size_ham(3)))

do i=1,size_ham(3)
   do j=1,Nspin
      allocate(anisotropy%E_line(j,i)%shell(indexNN(i)),anisotropy%mode_E_column(j,i)%shell(indexNN(i)))
   enddo
   call dissociate(anisotropy%E_line(:,i),indexNN(i),Nspin)
   call dissociate(anisotropy%mode_E_column(:,i),indexNN(i),Nspin)
enddo

end subroutine init_Energy_distrib

!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine that associates the Hamiltonians for each shells
!
!!!!!!!!!!!!!!!!!!!!!!!

subroutine associate_shell_ham_1D(ham_point,size_ham,ham_target,my_lattice,tableNN,indexNN)
implicit none
type(operator_real),intent(inout) :: ham_point(size_ham)
real(kind=8),target,intent(in) :: ham_target(:,:,:)
type(lattice), intent(in) :: my_lattice
integer, intent(in) :: tableNN(:,:,:,:,:,:) !!tableNN(4,N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),count(my_motif%i_mom)
integer, intent(in) :: indexNN(:)
integer, intent(in) :: size_ham
! internal variable
integer :: i,avant,n_atom_shell,Nspin,all_size(4)

all_size=shape(my_lattice%l_modes)
Nspin=product(all_size)

avant=0
do i=1,size_ham
!! i is the shell number
!! n_atom_shell is the number of atom in each shell

   n_atom_shell=indexNN(i)

   allocate(ham_point(i)%value(n_atom_shell,Nspin))
   allocate(ham_point(i)%line(n_atom_shell,Nspin))
   ham_point(i)%nline=Nspin
   ham_point(i)%ncolumn=Nspin
   ham_point(i)%line=0

   call dissociate(ham_point(i),Nspin,n_atom_shell)

   call associate_pointer(ham_point(i),ham_target(:,:,i),my_lattice,tableNN,indexNN(i),avant)

   avant=avant+indexNN(i)

enddo

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

size_line=shape(exchange%E_line)

do i=1,size_line(2)
   call associate_pointer(exchange%mode_E_column(:,i),spin,exchange%E_line(:,i),exchange%operators(i))
enddo

!!!!!!!!!!!!!!!
!
! Associate the DMI energy distribution
!
!!!!!!!!!!!!!!!

size_line=shape(DMI%E_line)

do i=1,size_line(2)
   call associate_pointer(DMI%mode_E_column(:,i),spin,DMI%E_line(:,i),DMI%operators(i))
enddo

!!!!!!!!!!!!!!!
!
! Associate the anisotropy energy distribution
!
!!!!!!!!!!!!!!!

size_line=shape(anisotropy%E_line)

do i=1,size_line(2)
   call associate_pointer(anisotropy%mode_E_column(:,i),spin,anisotropy%E_line(:,i),anisotropy%operators(i))
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
real(kind=8) :: E_ani,E_z,E_4,E_biq,E_dip,E_int
real(kind=8), allocatable :: E_DM(:),E_xch(:)
integer :: number_shell_exchange,number_shell_DMI
!   name of files
character(len=30) :: fname,rw_format

number_shell_exchange=size(exchange%operators)
number_shell_DMI=size(DMI%operators)

allocate(E_DM(number_shell_DMI),E_xch(number_shell_exchange))

write(rw_format,'( "(", I4, "(2x,f20.15))" )') number_shell_DMI+number_shell_exchange+5

E_DM=0.0d0
E_xch=0.0d0
E_ani=0.0d0
E_z=0.0d0
E_4=0.0d0
E_biq=0.0d0
E_dip=0.0d0
N=size(spin)

fname=convert('EnDistrib_',tag,'.dat')
io=open_file_write(fname)

do i=1,N

   E_xch=0.0d0
   do j=1,number_shell_exchange

   call local_energy_pointer(E_int,i,exchange%mode_E_column(i,j),exchange%E_line(i,j))
   E_xch(j)=E_xch(j)+E_int
   enddo

   E_DM=0.0d0
   do j=1,number_shell_DMI

   call local_energy_pointer(E_int,i,DMI%mode_E_column(i,j),DMI%E_line(i,j))
   E_DM(j)=E_DM(j)+E_int
   enddo

   call local_energy_pointer(E_int,i,anisotropy%mode_E_column(i,1),anisotropy%E_line(i,1))
   E_ani=E_int

   write(io,rw_format) E_DM,E_ani,E_z,E_4,E_biq,E_dip,E_xch
enddo
#ifdef CPP_BRUTDIP
!          if (i_dip) E_dip=dipole(i_x,i_y,i_z,i_m,spin,shape_spin,masque,my_lattice%boundary)
#else
!          if (i_dip) E_dip=fftdip(i_x,i_y,i_z,i_m)
#endif

!          Write(70,'(8(E20.12E3,2x))') E_DM,E_ani,E_z,E_4,E_biq,E_dip

Write(io,'(a)') '#1:E_DM   2:E_ani   3:E_z   4:E_4   5:E_biq   6:E_dip   7-N:E_xch   '

call close_file(fname,io)

end subroutine get_Energy_distrib

end module
