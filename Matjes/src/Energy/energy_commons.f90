module m_energy_commons
use m_basic_types, only : vec_point
use m_modes_variables, only : point_shell_mode
use m_derived_types, only : operator_real_order_N,point_shell_Operator,point_shell_Operator
use m_derived_types, only : lattice
use m_Hamiltonian_variables, only : coeff_ham_inter_spec,coeff_ham_inter_spec_pointer,shell_Ham,H_vois
use m_operator_pointer_utils
use m_lattice, only : my_order_parameters
use m_symmetry_operators
use m_total_Heisenberg_Ham
use m_summer_exp
!
! This has to be public so the variable can be accessed from elsewhere
! BUT it has to be protected so it can only be read and not written from outside the module
!
! coefficients of the Hamiltonian
type(coeff_ham_inter_spec_pointer), allocatable, public, protected, target :: Hamiltonians(:)
! total Hamiltonian
type(H_vois), public, protected, target, save :: total_hamiltonian
! total energy tensor
type(operator_real_order_N),public,protected,save :: energy
! translate the energy from a matrix * a coeff**N
real(kind=8), public, protected, allocatable, dimension(:,:,:) :: translator
! flags for the writting and reading of the Hamiltonian
logical :: write_ham=.false.
logical :: read_ham=.false.




private
public :: associate_energy_Hamiltonian,get_Hamiltonians,get_number_shell,get_total_Hamiltonian
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
use m_derived_types, only : simulation_parameters
use m_basic_types, only : site_Ham
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

! get the Heisenberg Hamiltonian
!call get_total_Heisenberg_Ham(fname,dim_ham,Ms)
!if (ham_tot_heisenberg%i_exist) number_hamiltonian=number_hamiltonian+1

! get the Hamiltonian for the summerfeld expansion
call get_Temperature_H(dim_ham)

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! put all the Hamiltonians into a big matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(Hamiltonians(number_hamiltonian))
number_hamiltonian=0

! Heisenberg
if (ham_tot_heisenberg%i_exist) then
  number_hamiltonian=number_hamiltonian+1
  Hamiltonians(number_hamiltonian)%name=ham_tot_heisenberg%name
  Hamiltonians(number_hamiltonian)%order=ham_tot_heisenberg%order
  Hamiltonians(number_hamiltonian)%onsite=.true.
  n_shell=size(ham_tot_heisenberg%ham)
  
  allocate(Hamiltonians(number_hamiltonian)%ham(n_shell))
  
  do i=1,n_shell
    Hamiltonians(number_hamiltonian)%ham(i)%Op_loc=>ham_tot_heisenberg%ham(i)%H
  enddo
  
endif

end subroutine get_Hamiltonians

!
! subroutines that prepare the total Hamiltonian
!
subroutine get_total_Hamiltonian(D,indexNN)
use m_Hamiltonian_variables, only : shell_Ham
use m_vector, only : norm
use m_io_utils
use m_convert
use m_io_files_utils
use m_setup_DM
implicit none
real(kind=8), intent(in) :: D(:,:,:)
integer, intent(in) :: indexNN(:,:)
! internal variables
integer :: n_ham,N_shell,dim_ham,n_ham_shell,N_atom_shell,n_order,i_order,avant
!slopes
integer :: i,j,l,k
! position in the matrices for the DMI
integer :: x_start,x_end,y_start,y_end
! order Ham
integer, allocatable :: order(:),shell_order_max(:,:)
integer :: max_order
! old stuff
integer :: io_input,n_DMI,N_voisin_DMI
type(shell_ham), allocatable, dimension(:) :: ham_DMI_Nshell_local
type(shell_ham), allocatable, dimension(:) :: ham_EM_Nshell_local
character(len=50) :: form

write(6,'(a)') 'Preparing the total Hamiltonian matrix'

n_ham=size(Hamiltonians)
dim_ham=size(Hamiltonians(1)%ham(1)%op_loc,1)
! number of Hamiltonian of order N
max_order=0
do i=1,n_ham
  if (Hamiltonians(i)%order.gt.max_order) max_order=Hamiltonians(i)%order
enddo
allocate(order(max_order))
order=0
do i=1,n_ham
  order( Hamiltonians(i)%order ) = order( Hamiltonians(i)%order ) + 1
enddo
n_order=count(order.ne.0)

write(6,'(a,I4,a)') 'found ',n_ham ,' Hamiltonian'
form=convert('(a,',max_order,'(2x,I3))')
write(6,form) 'order of the Hamitlonian',order
write(6,'(a)') ''

! set up the correct number of shells (always +1 for the onsite term always present in the total hamiltonians)
N_shell=size(indexNN,1)
allocate(total_hamiltonian%shell_num(N_shell+1),total_hamiltonian%num(N_shell+1))
! always the onsite term
total_hamiltonian%num=1
allocate(total_hamiltonian%shell_num(1)%order(1),total_hamiltonian%shell_num(1)%num(1))
total_hamiltonian%shell_num(1)%num=2
! only one atom the site 0
allocate(total_hamiltonian%shell_num(1)%order(1)%atom(1))

! allocate Hamiltonian of the anisotropy and the Zeeman
allocate(total_hamiltonian%shell_num(1)%order(1)%atom(1)%H(dim_ham,dim_ham))
total_hamiltonian%shell_num(1)%order(1)%line=dim_ham
total_hamiltonian%shell_num(1)%order(1)%column=dim_ham
total_hamiltonian%shell_num(1)%order(1)%atom(1)%H=0.0d0

! load the number of neighbors in eachs shell in the total_Hamiltonians
do i=1,N_shell
  total_hamiltonian%num(i+1)=indexNN(i,1)
enddo

write(6,'(a,I4,a)') 'Total Hamiltonian has ',N_shell+1 ,' Shells'
form=convert('(a,',N_shell+1,'(2x,I3))')
write(6,form) 'Number of sites per shell:',total_hamiltonian%num
write(6,'(a)') ''

!
! find the shells with order N
!

!get shell_order_max
!  shell_order_max(:,1) maximum Hamiltonian order for each shell
!  shell_order_max(:,2) number of Hamiltonian per shell
!  Hamiltonians(i)%onsite=.True. means that the shell 0 is included in the Hamiltonian
allocate(shell_order_max(N_shell+1,2))
shell_order_max=0
do i=1,n_ham
  n_ham_shell=size(Hamiltonians(i)%ham)
  i_order=Hamiltonians(i)%order
  do j=1,n_ham_shell
    if (Hamiltonians(i)%onsite) then
      if (shell_order_max(j,1).lt.i_order) shell_order_max(j,1)=i_order
      shell_order_max(j,2)=shell_order_max(j,2)+1
    else
      if (shell_order_max(j+1,1).lt.i_order) shell_order_max(j+1,1)=i_order
      shell_order_max(j+1,2)=shell_order_max(j+1,2)+1
    endif
  enddo
enddo
form=convert('(a,',N_shell+1,'I4)')
write(6,form) 'Maximum Hamiltonianss order per shell', shell_order_max(:,1)
write(6,form) 'number of Hamiltonians per shell', shell_order_max(:,2)

!allocate correct sizes for total_hamiltonian
do i=2,N_shell+1
  n_order=shell_order_max(i,1)
  ! if you have a Hamiltonian of order 3 on one shell, you must an Hamiltonian of order 2
  allocate(total_hamiltonian%shell_num(i)%order(n_order-1),total_hamiltonian%shell_num(i)%num(n_order-1))
  total_hamiltonian%shell_num(i)%num=n_order-1

  N_atom_shell=total_hamiltonian%num(i)
  do i_order=2,n_order
    allocate(total_hamiltonian%shell_num(i)%order(i_order-1)%atom( N_atom_shell ))
    total_hamiltonian%shell_num(i)%order(i_order-1)%line=dim_ham**(i_order-1)
    total_hamiltonian%shell_num(i)%order(i_order-1)%column=dim_ham
    do j=1,total_hamiltonian%num(i)
      allocate(total_hamiltonian%shell_num(i)%order(i_order-1)%atom(j)%H( dim_ham , dim_ham**( i_order-1) ))
      total_hamiltonian%shell_num(i)%order(i_order-1)%atom(j)%H=0.0d0
    enddo
  enddo
enddo

!
! I do not know the order of interaction of each shell. Therefore the shell Hamiltonian must be specified first
!

!
! prepare all the onsite Hamiltonians
!

do i=1,n_ham

!
! Heisenberg part
!

  if ('heisenberg'.eq.trim(Hamiltonians(i)%name)) then
    n_DMI=get_number_DMI('input')
    N_shell=size(Hamiltonians(i)%ham)

    ! ham_DMI_Nshell_local( Hamiltonian column , Hamiltonian row , number of atom in shell i , shell i )

    allocate(ham_DMI_Nshell_local(N_shell-1))

    do j=2,N_shell
       N_voisin_DMI=indexNN(j-1,1)
       allocate(ham_DMI_Nshell_local(j-1)%atom(N_voisin_DMI))
       ham_DMI_Nshell_local(j-1)%line=dim_ham
       ham_DMI_Nshell_local(j-1)%column=dim_ham
       do k=1,size(ham_DMI_Nshell_local(j-1)%atom)
          allocate(ham_DMI_Nshell_local(j-1)%atom(k)%H(dim_ham,dim_ham))
          ham_DMI_Nshell_local(j-1)%atom(k)%H=0.0d0
       enddo
    enddo

    call get_borders('magnetic',x_start,x_end,'magnetic',y_start,y_end,my_order_parameters)

    ! Hamiltonians(i)%ham is just the shell number. The number of atom in each shell is not taken into account here
    do j=2,size(Hamiltonians(i)%ham)
      if (j.le.n_DMI+1) then
         do k=1,size(ham_DMI_Nshell_local(j-1)%atom)
            call convoluate_Op_SOC_vector(D(k,:,1),Hamiltonians(i)%ham(j)%op_loc(x_start:x_end,y_start:y_end),ham_DMI_Nshell_local(j-1)%atom(k)%H(x_start:x_end,y_start:y_end))
         enddo
      else
         do k=1,size(ham_DMI_Nshell_local(j-1)%atom)
            write(*,*) j,k
            ham_DMI_Nshell_local(j-1)%atom(k)%H=Hamiltonians(i)%ham(j)%op_loc
         enddo
      endif

    enddo
    ! The number of atom in the shell is then in ham_DMI_Nshell_local

  endif
!
! Magnetoelectric part
!
  if ('magnetoelectric'.eq.trim(Hamiltonians(i)%name)) then
    N_shell=size(Hamiltonians(i)%ham)
    n_ham_shell=sum(indexNN(1:N_shell,1))

    allocate(ham_EM_Nshell_local(N_shell))

    do j=1,size(ham_EM_Nshell_local)
       N_voisin_DMI=indexNN(j,1)
       allocate(ham_EM_Nshell_local(j)%atom(N_voisin_DMI))
       ham_EM_Nshell_local(j)%line=dim_ham**2
       ham_EM_Nshell_local(j)%column=dim_ham
       do k=1,size(ham_EM_Nshell_local(j)%atom)
          allocate(ham_EM_Nshell_local(j)%atom(k)%H(dim_ham,dim_ham**2))
          ham_EM_Nshell_local(j)%atom(k)%H=0.0d0
       enddo
    enddo

    call get_borders('magnetic',x_start,x_end,'Efield',y_start,y_end,my_order_parameters)

    ! Hamiltonians(i)%ham is just the shell number. The number of atom in each shell is not taken into account here
    do j=1,size(Hamiltonians(i)%ham)
      do k=1,size(ham_EM_Nshell_local(j)%atom)
        call convoluate_Op_sym(j,k,Hamiltonians(i)%ham(j)%op_loc,ham_EM_Nshell_local(j)%atom(k)%H,x_start,x_end,y_start,y_end)
      enddo
    enddo

  endif

enddo

! put all the different parts of the Hamiltonian in the big one
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! first the onsite terms
do i=1,n_ham

  if ('heisenberg'.eq.trim(Hamiltonians(i)%name)) then
    total_hamiltonian%shell_num(1)%order(1)%atom(1)%H=total_hamiltonian%shell_num(1)%order(1)%atom(1)%H+Hamiltonians(i)%ham(1)%op_loc
    do j=2,size(Hamiltonians(i)%ham)
      do l=1,size(ham_DMI_Nshell_local(j-1)%atom)
        total_hamiltonian%shell_num(j)%order(1)%atom(l)%H=total_hamiltonian%shell_num(j)%order(1)%atom(l)%H+ham_DMI_Nshell_local(j-1)%atom(l)%H
      enddo
    enddo
  endif

  if ('magnetoelectric'.eq.trim(Hamiltonians(i)%name)) then
    do j=1,size(ham_EM_Nshell_local) ! carefull the shell 1 is for the anisotropy (I count the shell in the EM Hamiltonian here)
      total_hamiltonian%shell_num(j+1)%num(2)=3
      do l=1,total_hamiltonian%num(j+1)   ! so one needs to look at j+1 since j=1 is the anisotropy tensor of order 2
        total_hamiltonian%shell_num(j+1)%order(2)%atom(l)%H=total_hamiltonian%shell_num(j+1)%order(2)%atom(l)%H+ham_EM_Nshell_local(j)%atom(l)%H
      enddo
    enddo
  endif
  
enddo

! prepare the translator for the summerfeld expansion
if (temperature_strasbourg%i_exist) then
    N_shell=size(total_hamiltonian%shell_num)
    n_ham_shell=sum(indexNN(1:N_shell-1,1))+1
    allocate(translator(dim_ham,dim_ham,n_ham_shell))
    translator=0.0d0
    translator(:,:,1)=temperature_strasbourg%ham(1)%H
    avant=1
    do i=2,N_shell
      do j=1,indexNN(i-1,1)
        translator(:,:,avant+j)=temperature_strasbourg%ham(i)%H
      enddo
      avant=avant+indexNN(i-1,1)
    enddo
endif

form=convert('(',dim_ham,'(2x,f12.8))')
io_input=open_file_read('input')
call get_parameter(io_input,'input','write_ham',write_ham)
call close_file('input',io_input)
if (write_ham) call rw_hamiltonian('Hamiltonian.dat')

write(6,'(a)') ''
write(6,'(a)') 'Total Hamiltonian'
do i=1,size(total_hamiltonian%shell_num)
  write(6,'(a,2x,I3)') 'Hamiltonian for shell ',i

  do j=1,size(total_hamiltonian%shell_num(i)%order)
    write(6,'(a,2x,I3)') 'Hamiltonian of order ',j

    do k=1,total_hamiltonian%num(i)

        write(6,'(a,2x,I3)') 'Bond', k

      do l=1,total_hamiltonian%shell_num(i)%order(j)%line
        write(6,form) total_hamiltonian%shell_num(i)%order(j)%atom(k)%H(:,l)
      enddo

    enddo
    write(6,'(a)') ''

  enddo
  write(6,'(a)') ''
enddo

end subroutine get_total_Hamiltonian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! count the number of shell needed for the simulation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer function get_number_shell()
implicit none
! internal
integer :: N_coeff
integer :: i

get_number_shell=0

do i=1,size(Hamiltonians)
! obtain the size of all matrices
   N_coeff=size(Hamiltonians(i)%ham)
   if (N_coeff.gt.get_number_shell) get_number_shell=N_coeff
enddo

if (get_number_shell.eq.0) then
   write(6,'(a)') 'WARNING!!!! No shells were found in the Hamiltonian (H=0)'
else
   ! the onsite term is always preset in the Hamiltonian so we take it out for the table of distance and the table of neighbor
   get_number_shell=get_number_shell-1
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
integer :: i,j,k,column,l
integer :: io
character(len=50) :: form

io=open_file_write(fname)

do i=1,size(total_hamiltonian%shell_num)
  write(io,'(a,2x,I3)') 'Hamiltonian for shell ',i
  do j=1,size(total_hamiltonian%shell_num(i)%order)
    write(io,'(a,2x,I3)') 'Hamiltonian of order ',j

    do k=1,total_hamiltonian%num(i)

      write(io,'(a,2x,I3)') 'Bond', k

      column=total_hamiltonian%shell_num(i)%order(j)%column
      form=convert('(',column,'(E12.6,2x))')

      do l=1,total_hamiltonian%shell_num(i)%order(j)%line
        write(io,form) total_hamiltonian%shell_num(i)%order(j)%atom(k)%H(:,l)
      enddo

      write(io,'(a)') ''
    enddo
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

all_size=shape(my_lattice%ordpar%l_modes)
Nspin=product(all_size)
shape_tableNN=shape(tableNN)

allocate(energy%value(shape_tableNN(2)+1,Nspin))
allocate(energy%line(shape_tableNN(2)+1,Nspin))
energy%nline=Nspin
energy%ncolumn=Nspin
energy%line=0

call associate_pointer(energy,total_hamiltonian,my_lattice,tableNN,indexNN)

end subroutine associate_energy_Hamiltonian

end module m_energy_commons
