module m_energyfield
use m_derived_types
use m_operator_pointer_utils

! Exchange energy tensor
type(operator_real),allocatable :: exchange(:)
! DMI energy tensor
type(operator_real),allocatable :: DMI(:)
! anisotropy energy tensor
type(operator_real) :: anisotropy
! Zeeman energy tensor
type(operator_real) :: Zeeman

type(point_shell_Operator), allocatable, dimension(:) :: E_line
type(point_shell_mode), allocatable, dimension(:) :: mode_E_column

private
public :: init_Energy_distrib,get_Energy_distrib
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
integer :: size_ham(3),Nspin,all_size(4),shape_tableNN(6)
! slope variables
integer :: i

all_size=size(my_lattice%l_modes)
size_ham=size(Hamiltonian%exchange)
Nspin=product(all_size)
shape_tableNN=shape(tableNN)
allocate(exchange(size_ham(3)))

write(*,*) Nspin
pause
do i=1,size_ham(3)
!! here is the shell number
   allocate(exchange(i)%value(Nspin,Nspin))
   exchange(i)%nline=Nspin
   exchange(i)%ncolumn=Nspin

   call dissociate(exchange(i))

   call associate_pointer(exchange(i),Hamiltonian%exchange(:,:,i),my_lattice,tableNN,indexNN(i))

enddo

size_ham=size(Hamiltonian%DMI)
allocate(DMI(size_ham(3)))

stop

end subroutine init_Energy_distrib

!!!!!!!!!!!!!!!!!!!!!!!
!
! Print the energy density distribution interaction and shell resolve into a file
!
!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_Energy_distrib(tag,spin)
use m_io_files_utils
use m_convert
implicit none
! input
type(vec_point),intent(in) :: spin(:)
integer, intent(in) :: tag
! internals
integer :: i,N,io
real(kind=8) :: E_ani,E_z,E_4,E_biq,E_dip
real(kind=8), allocatable :: E_DM(:),E_xch(:)
!   name of files
character(len=30) :: fname

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

Write(io,'(a)') '#1:E_DM   2:E_ani   3:E_z   4:E_4   5:E_biq   6:E_dip   7-N:E_xch   '

do i=1,N

#ifdef CPP_BRUTDIP
!          if (i_dip) E_dip=dipole(i_x,i_y,i_z,i_m,spin,shape_spin,masque,my_lattice%boundary)
#else
!          if (i_dip) E_dip=fftdip(i_x,i_y,i_z,i_m)
#endif

!          Write(70,'(8(E20.12E3,2x))') E_DM,E_ani,E_z,E_4,E_biq,E_dip

enddo

call close_file(fname,io)

end subroutine get_Energy_distrib

end module
