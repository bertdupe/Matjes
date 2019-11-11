module m_energy_commons
use m_derived_types, only : coeff_ham_inter_spec,coeff_ham_inter_spec_pointer,operator_real,Coeff_Ham,shell_Ham,point_shell_Operator,vec_point,point_shell_mode,lattice
use m_operator_pointer_utils
use m_lattice, only : my_order_parameters
use m_couplage_ME
use m_exchange_heisenberg
use m_anisotropy_heisenberg
use m_zeeman
use m_symmetry_operators
!
! This has to be public so the variable can be accessed from elsewhere
! BUT it has to be protected so it can only be read and not written from outside the module
!
! coefficients of the Hamiltonian
type(coeff_ham_inter_spec_pointer), allocatable, public, protected, target :: Hamiltonians(:)
!type(Coeff_Ham),public,protected,target,save :: Hamiltonian
! total Hamiltonian
type(shell_Ham), public, protected, target, save, allocatable, dimension(:) :: total_hamiltonian
! total energy tensor
type(operator_real),public,protected,save :: energy
! flags for the writting and reading of the Hamiltonian
logical :: write_ham=.false.
logical :: read_ham=.false.




private
public :: associate_energy_Hamiltonian,get_Hamiltonians,get_number_shell,get_total_Hamiltonian,get_E_line
contains

!!!!!!!!!!!!!!!!!!!!!!!
!  Routine that prepares the LOCAL Hamiltonian for a particular geometry
! The form of the Hamiltonian is
!
!  ( H_00  ...           )  (Site_0)
!  ( H_10  H_11          )  (Site_1)
!  (  .          .       )  (Site_2)
!  (  .             .    )  (...   )
!  (  .               .  )  (...   )
!
!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_Hamiltonians(fname,Ms,dim_Ham)
use m_derived_types, only : site_Ham
use m_io_files_utils
use m_io_utils
use m_constants, only : mu_B
implicit none
character(len=*), intent(in) :: fname
real(kind=8), intent(in) :: Ms
integer, intent(in) :: dim_Ham
! internal
integer :: number_hamiltonian,n_shell,i

number_hamiltonian=0

! get the anisotropy Hamiltonian
call get_ham_anisotropy(fname,dim_ham)
if (anisotropy%i_exist) number_hamiltonian=number_hamiltonian+1

! get the zeeman
call get_ham_zeeman(fname,dim_ham,Ms)
if (zeeman%i_exist) number_hamiltonian=number_hamiltonian+1

! get the exchange Hamiltonian
call get_ham_exchange(fname,dim_ham)
if (exchange%i_exist) number_hamiltonian=number_hamiltonian+1

! get the magnetoelectric coefficients
call get_ham_ME(fname,dim_ham)
if (ME%i_exist) number_hamiltonian=number_hamiltonian+1

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! put all the Hamiltonians into a big matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(Hamiltonians(number_hamiltonian))
number_hamiltonian=0
! the anisotropy first
if (anisotropy%i_exist) then
  number_hamiltonian=number_hamiltonian+1
  Hamiltonians(number_hamiltonian)%name=anisotropy%name
  n_shell=size(anisotropy%ham)
  allocate(Hamiltonians(number_hamiltonian)%ham(n_shell))
  do i=1,n_shell
     Hamiltonians(number_hamiltonian)%ham(i)%Op_loc=>anisotropy%ham(i)%H
  enddo
endif

! the exchange
if (exchange%i_exist) then
  number_hamiltonian=number_hamiltonian+1
  Hamiltonians(number_hamiltonian)%name=exchange%name
  n_shell=size(exchange%ham)
  allocate(Hamiltonians(number_hamiltonian)%ham(n_shell))
  do i=1,n_shell
     Hamiltonians(number_hamiltonian)%ham(i)%Op_loc=>exchange%ham(i)%H
  enddo
endif

! The magnetoelectric coupling
if (ME%i_exist) then
  number_hamiltonian=number_hamiltonian+1
  Hamiltonians(number_hamiltonian)%name=ME%name
  n_shell=size(ME%ham)
  allocate(Hamiltonians(number_hamiltonian)%ham(n_shell))
  do i=1,n_shell
     Hamiltonians(number_hamiltonian)%ham(i)%Op_loc=>ME%ham(i)%H
  enddo
endif

! The Zeeman
if (Zeeman%i_exist) then
  number_hamiltonian=number_hamiltonian+1
  Hamiltonians(number_hamiltonian)%name=zeeman%name
  allocate(Hamiltonians(number_hamiltonian)%ham(1))
  Hamiltonians(number_hamiltonian)%ham(1)%Op_loc=>zeeman%ham(1)%H
endif

end subroutine get_Hamiltonians

!
! subroutines that prepare the total Hamiltonian
!
subroutine get_total_Hamiltonian(D,indexNN)
use m_derived_types, only : shell_Ham
use m_vector, only : norm
use m_io_utils
use m_io_files_utils
use m_setup_DM
implicit none
real(kind=8), intent(in) :: D(:,:,:)
integer, intent(in) :: indexNN(:,:)
! internal variables
integer :: n_ham,N_shell,dim_ham,n_ham_shell,N_atom_shell
integer, allocatable :: size_hams(:)
!slopes
integer :: i,j,l,k
! position in the matrices for the DMI
integer :: x_start,x_end,y_start,y_end
! old stuff
integer :: N_tot_vec_shell,N_max_shell,N_max_nei,io_input,n_DMI,N_voisin_DMI
integer :: dim_D(3)
type(shell_ham), allocatable, dimension(:) :: ham_DMI_Nshell_local
type(shell_ham), allocatable, dimension(:) :: ham_EM_Nshell_local

! obtain the size of all matrices
n_ham=size(Hamiltonians)
allocate(size_hams(n_ham))
size_hams=0
do i=1,n_ham
   size_hams(i)=size(Hamiltonians(i)%ham)
enddo
N_max_shell=maxval(size_hams)
dim_D=shape(D)
! the number of coefficicient must match the number of DMI shells
N_tot_vec_shell=size(indexNN(:,1),1)
N_max_nei=maxval(indexNN(:,1),1)

do i=1,n_ham

!
! Exchange and DMI part part
!

  if ('exchange'.eq.trim(Hamiltonians(i)%name)) then
    n_DMI=get_number_DMI('input')
    N_shell=size(Hamiltonians(i)%ham)
    n_ham_shell=sum(indexNN(1:N_shell,1))
    dim_ham=size(Hamiltonians(i)%ham(1)%op_loc,1)

    ! ham_DMI_Nshell_local( Hamiltonian column , Hamiltonian row , number of atom in shell i , shell i )

    allocate(ham_DMI_Nshell_local(N_shell))

    do j=1,size(ham_DMI_Nshell_local)
       N_voisin_DMI=indexNN(j,1)
       allocate(ham_DMI_Nshell_local(j)%atom(N_voisin_DMI))
       do k=1,size(ham_DMI_Nshell_local(j)%atom)
          allocate(ham_DMI_Nshell_local(j)%atom(k)%H(dim_ham,dim_ham))
          ham_DMI_Nshell_local(j)%atom(k)%H=0.0d0
       enddo
    enddo

    call get_borders('magnetic',x_start,x_end,'magnetic',y_start,y_end,my_order_parameters)

    ! Hamiltonians(i)%ham is just the shell number. The number of atom in each shell is not taken into account here
    do j=1,size(Hamiltonians(i)%ham)
      if (j.le.n_DMI) then
         do k=1,size(ham_DMI_Nshell_local(j)%atom)
            call convoluate_Op_SOC_vector(D(k,:,1),Hamiltonians(i)%ham(j)%op_loc(x_start:x_end,y_start:y_end),ham_DMI_Nshell_local(j)%atom(k)%H(x_start:x_end,y_start:y_end))
         enddo
      else
         do k=1,size(ham_DMI_Nshell_local(j)%atom)
            ham_DMI_Nshell_local(j)%atom(k)%H=Hamiltonians(i)%ham(j)%op_loc
         enddo
      endif

    enddo
    ! The number of atom in the shell is then in ham_DMI_Nshell_local

  endif

!
! Magnetoelectric part
!
  if ('magnetoelectric'.eq.trim(Hamiltonians(i)%name)) then
    n_DMI=get_number_EM_DMI('input')
    N_shell=size(Hamiltonians(i)%ham)
    n_ham_shell=sum(indexNN(1:N_shell,1))
    dim_ham=size(Hamiltonians(i)%ham(1)%op_loc,1)

    allocate(ham_EM_Nshell_local(N_shell))

    do j=1,size(ham_EM_Nshell_local)
       N_voisin_DMI=indexNN(j,1)
       allocate(ham_EM_Nshell_local(j)%atom(N_voisin_DMI))
       do k=1,size(ham_EM_Nshell_local(j)%atom)
          allocate(ham_EM_Nshell_local(j)%atom(k)%H(dim_ham,dim_ham))
          ham_EM_Nshell_local(j)%atom(k)%H=0.0d0
       enddo
    enddo

    call get_borders('magnetic',x_start,x_end,'Efield',y_start,y_end,my_order_parameters)

    ! Hamiltonians(i)%ham is just the shell number. The number of atom in each shell is not taken into account here
    do j=1,size(Hamiltonians(i)%ham)
      if (j.le.n_DMI) then
         do k=1,size(ham_EM_Nshell_local(j)%atom)
            call convoluate_Op_SOC_vector(D(k,:,1),Hamiltonians(i)%ham(j)%op_loc(x_start:x_end,y_start:y_end),ham_EM_Nshell_local(j)%atom(k)%H(x_start:x_end,y_start:y_end))
            call convoluate_Op_SOC_vector(D(k,:,1),Hamiltonians(i)%ham(j)%op_loc(y_start:y_end,x_start:x_end),ham_EM_Nshell_local(j)%atom(k)%H(y_start:y_end,x_start:x_end))
         enddo
      else
         do k=1,size(ham_DMI_Nshell_local(j)%atom)
            ham_EM_Nshell_local(j)%atom(k)%H=Hamiltonians(i)%ham(j)%op_loc
         enddo
      endif
    enddo

  endif

enddo


allocate(total_hamiltonian(N_tot_vec_shell+1))
allocate(total_hamiltonian(1)%atom(1))
allocate(total_hamiltonian(1)%atom(1)%H(dim_ham,dim_ham))
total_hamiltonian(1)%atom(1)%H=0.0d0

do i=1,N_tot_vec_shell
   N_atom_shell=indexNN(i,1)
   allocate(total_hamiltonian(i+1)%atom(N_atom_shell))
   do j=1,N_atom_shell
      allocate(total_hamiltonian(i+1)%atom(j)%H(dim_ham,dim_ham))
      total_hamiltonian(i+1)%atom(j)%H=0.0d0
   enddo
enddo

! put all the different parts of the Hamiltonian in the big one
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! first the onsite terms

do i=1,n_ham

  if ('zeeman'.eq.trim(Hamiltonians(i)%name)) call add(total_hamiltonian(1)%atom(1)%H,Hamiltonians(i)%ham(1)%Op_loc)
  if ('anisotropy'.eq.trim(Hamiltonians(i)%name)) call add(total_hamiltonian(1)%atom(1)%H,Hamiltonians(i)%ham(1)%Op_loc)

  if ('exchange'.eq.trim(Hamiltonians(i)%name)) then
    do j=1,size(ham_DMI_Nshell_local)
      do l=1,size(total_hamiltonian(j+1)%atom)
        total_hamiltonian(j+1)%atom(l)%H=total_hamiltonian(j+1)%atom(l)%H+ham_DMI_Nshell_local(j)%atom(l)%H
      enddo
    enddo
  endif

  if ('magnetoelectric'.eq.trim(Hamiltonians(i)%name)) then
    do j=1,size(ham_EM_Nshell_local)
      do l=1,size(total_hamiltonian(j+1)%atom)
        total_hamiltonian(j+1)%atom(l)%H=total_hamiltonian(j+1)%atom(l)%H+ham_EM_Nshell_local(j)%atom(l)%H
      enddo
    enddo
  endif
enddo

io_input=open_file_read('input')
call get_parameter(io_input,'input','write_ham',write_ham)
call close_file('input',io_input)
if (write_ham) call rw_hamiltonian('Hamiltonian.dat')

#ifdef CPP_DEBUG
do i=1,N_max_shell+1
   write(*,*) 'shell', i

   do j=1,size(total_hamiltonian(i)%atom)
      write(6,*) (total_hamiltonian(i)%atom(j)%H(:,k),k=1,dim_ham)
   enddo

!do i=1,size(Hamiltonians)
!   do j=1,size(Hamiltonians(i)%ham)
!      write(6,*) (Hamiltonians(i)%ham(j)%Op_loc(:,k),k=1,dim_ham)
!   enddo

enddo
stop 'toto'
#endif

end subroutine get_total_Hamiltonian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! count the number of shell needed for the simulation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer function get_number_shell()
implicit none
! internal
integer :: N_coeff
integer :: i

get_number_shell=size(Hamiltonians(1)%ham)

do i=1,size(Hamiltonians)
! obtain the size of all matrices
   N_coeff=size(Hamiltonians(i)%ham)
   if (N_coeff.gt.get_number_shell) get_number_shell=N_coeff
enddo

if (get_number_shell.eq.0) then
   write(6,'(a)') 'WARNING!!!! No shells were found in the Hamiltonian (H=0)'
else
   write(6,'(I6,a/)') N_coeff,' shells were found in the Hamiltonian'
endif

end function get_number_shell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! add a pointer of a good size to the matrix of the total hamiltonian
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine add(H,pointer)
implicit none
real(kind=8), intent(inout) :: H(:,:)
real(kind=8), intent(in) :: pointer(:,:)
!internal

if (size(H).ne.size(pointer)) stop 'ERROR in energy_commons ADD function'

H=H+pointer

end subroutine add

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! smal subroutine that reads and writes the total Hamiltonian
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rw_hamiltonian(fname)
use m_io_files_utils
use m_convert
implicit none
character(len=*), intent(in) :: fname
! internal
integer :: i,j,k,dim_ham
integer :: io
character(len=30) :: format

io=open_file_write(fname)

do i=1,size(total_hamiltonian)
   write(io,'(a,2x,I4)') 'shell',i
   do j=1,size(total_hamiltonian(i)%atom)
      dim_ham=size(total_hamiltonian(i)%atom(j)%H,2)
      format=convert('(',dim_ham,'(E12.6,2x))')

      write(io,format) (total_hamiltonian(i)%atom(j)%H(:,k),k=1,dim_ham)

      write(io,'(a)') ''
   enddo
   write(io,'(a)') ''
enddo

call close_file(fname,io)


end subroutine rw_hamiltonian










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! prepare the Hamiltonian matrix and the pointers for the matrix Hamiltonian static association
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine associate_energy_Hamiltonian(my_lattice,tableNN,indexNN)
implicit none
type(lattice), intent(in) :: my_lattice
integer, intent(in) :: tableNN(:,:,:,:,:,:) !!tableNN(4,N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),count(my_motif%i_mom)
integer, intent(in) :: indexNN(:)
! internal
! slope of the sums
integer :: Nspin,all_size(4),shape_tableNN(6)

all_size=shape(my_lattice%l_modes)
Nspin=product(all_size)
shape_tableNN=shape(tableNN)

allocate(energy%value(shape_tableNN(2)+1,Nspin))
allocate(energy%line(shape_tableNN(2)+1,Nspin))
energy%nline=Nspin
energy%ncolumn=Nspin
energy%line=0

call dissociate(energy%value,shape_tableNN(2)+1,Nspin)

call associate_pointer(energy,total_hamiltonian,my_lattice,tableNN,indexNN)

! old version
!!!call associate_pointer(energy,Hamiltonian%total_shell,my_lattice,tableNN,indexNN)

end subroutine associate_energy_Hamiltonian

subroutine get_E_line(E_line,mode_E_column,spin)
implicit none
type(point_shell_Operator), intent(inout) :: E_line(:)
type(point_shell_mode),intent(inout) :: mode_E_column(:)
type(vec_point),target,intent(in) :: spin(:)
! internal variables
integer :: shape_energy(2)
! slope variables
integer :: i

shape_energy=shape(energy%value)

do i=1,shape_energy(2)
   allocate(E_line(i)%shell(shape_energy(1)))
   allocate(mode_E_column(i)%shell(shape_energy(1)))
enddo
call dissociate(E_line,shape_energy(1),shape_energy(2))
call dissociate(mode_E_column,shape_energy(1),shape_energy(2))


call associate_pointer(mode_E_column,spin,E_line,energy)

end subroutine get_E_line

end module m_energy_commons
