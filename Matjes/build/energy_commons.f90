module m_energy_commons
use m_derived_types, only : operator_real,Coeff_Ham,shell_Ham,point_shell_Operator,vec_point,point_shell_mode
use m_operator_pointer_utils
! coefficients of the Hamiltonian
type(Coeff_Ham),public,protected,target,save :: Hamiltonian
! total energy tensor
type(operator_real),public,protected,save :: energy
! Exchange energy tensor
type(operator_real),public,protected,save :: exchange
! DMI energy tensor
type(operator_real),public,protected,save :: DMI
! anisotropy energy tensor
type(operator_real),public,protected,save :: anisotropy
! Zeeman energy tensor
type(operator_real),public,protected,save :: Zeeman

interface number_nonzero_coeff
  module procedure number_nonzero_coeff_1d
  module procedure number_nonzero_coeff_2d
end interface number_nonzero_coeff

private
public :: associate_energy_Hamiltonian,get_Hamiltonians,get_number_shell,get_total_Hamiltonian,associate_line_Hamiltonian
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

subroutine get_Hamiltonians(fname,B)
use m_derived_types, only : site_Ham
use m_io_files_utils
use m_io_utils
use m_constants, only : mu_B
implicit none
character(len=*), intent(in) :: fname
real(kind=8), intent(in) :: B(:)
! internal
integer :: io_param,neighbor_Jij,neighbor_DM,neighbor_ani
! maximum number of shell to consider
integer :: N_coeff_Exch,N_coeff_DMI,N_coeff_ani
! exchange and DMI
real(kind=8), allocatable :: J_ij(:),D_ij(:),ani(:,:)
! slope
integer :: i,j

neighbor_Jij=0
neighbor_DM=0
neighbor_ani=0
io_param=open_file_read(fname)

! count the neighbors for the exchange
neighbor_Jij=count_variables(io_param,'J_','input')
allocate(J_ij(neighbor_Jij))
J_ij=0.0d0

neighbor_DM=count_variables(io_param,'DM_','input')
allocate(D_ij(neighbor_DM))
D_ij=0.0d0

neighbor_ani=count_variables(io_param,'ani_','input')
if (neighbor_ani.eq.0) neighbor_ani=1
allocate(ani(3,neighbor_ani))
ani=0.0d0

! get the exchange, DMI, anisotropy
call get_coeff(io_param,'input','J_',J_ij)
N_coeff_Exch=number_nonzero_coeff(J_ij,'exchange')
if (N_coeff_Exch.gt.0) Hamiltonian%i_exch=.true.

call get_coeff(io_param,'input','DM_',D_ij)
N_coeff_DMI=number_nonzero_coeff(D_ij,'DMI')
if (N_coeff_DMI.gt.0) Hamiltonian%i_DM=.true.

call get_parameter(io_param,'input','ani_',3,ani(:,1))
N_coeff_ani=number_nonzero_coeff(ani,'anisotropy')
if (N_coeff_ani.gt.0) Hamiltonian%i_ani=.true.
if (N_coeff_ani.eq.0) N_coeff_ani=1

call close_file(fname,io_param)

! allocate the variables
allocate(Hamiltonian%exchange(3,3,N_coeff_Exch))
allocate(Hamiltonian%ani(3,3,N_coeff_ani))
allocate(Hamiltonian%DMI(3,3,N_coeff_DMI))
allocate(Hamiltonian%Zeeman(3,3))

! initialize variables
Hamiltonian%exchange=0.0d0
Hamiltonian%ani=0.0d0
Hamiltonian%DMI=0.0d0
Hamiltonian%Zeeman=0.0d0

! set up the anisotropy

do i=1,N_coeff_ani
   do j=1,3
      Hamiltonian%ani(j,j,i)=ani(j,i)*Hamiltonian%c_ani
   enddo
enddo

! set up the Zeeman energy
do j=1,3
   Hamiltonian%Zeeman(j,j)=B(j)
enddo

! set up the magnetic exchange interaction
do i=1,N_coeff_Exch
   do j=1,3
      Hamiltonian%exchange(j,j,i)=J_ij(i)*Hamiltonian%c_Ji
   enddo
enddo

! set up the DMI
do i=1,N_coeff_DMI
   Hamiltonian%DMI(2,1,i)=-D_ij(i)*Hamiltonian%c_DM
   Hamiltonian%DMI(3,1,i)=D_ij(i)*Hamiltonian%c_DM
   Hamiltonian%DMI(3,2,i)=-D_ij(i)*Hamiltonian%c_DM

   Hamiltonian%DMI(1,2,i)=-Hamiltonian%DMI(2,1,i)
   Hamiltonian%DMI(2,3,i)=-Hamiltonian%DMI(3,2,i)
   Hamiltonian%DMI(1,3,i)=-Hamiltonian%DMI(3,1,i)
enddo

deallocate(J_ij,D_ij,ani)

end subroutine get_Hamiltonians

!
! subroutines that prepare the total Hamiltonian
!
subroutine get_total_Hamiltonian(D,indexNN)
use m_derived_types, only : shell_Ham
implicit none
real(kind=8), intent(in) :: D(:,:,:)
integer, intent(in) :: indexNN(:,:)
! internal variables
integer :: N_coeff_ani,N_coeff_DMI,N_coeff_Exch,N_Dim_Htot
integer :: N_tot_vec_shell,N_max_shell,N_max_nei
integer :: i,j,l

! obtain the size of all matrices
N_coeff_Exch=size(Hamiltonian%exchange,3)
N_coeff_ani=size(Hamiltonian%ani,3)
if (N_coeff_ani.eq.0) N_coeff_ani=1
N_coeff_DMI=size(Hamiltonian%DMI,3)
N_max_shell=max(N_coeff_Exch,N_coeff_ani,N_coeff_DMI)

! the number of coefficicient must match the number of DMI shells
N_tot_vec_shell=size(indexNN(:,1),1)
N_max_nei=maxval(indexNN(:,1),1)

! total number of matrices to take into account into the total energy
N_Dim_Htot=1    ! that is to take into account the anisotropy
do i=1,N_tot_vec_shell

   if (i.le.N_coeff_DMI) then
      if ( all( abs(Hamiltonian%DMI(:,:,i)).lt.1.0d-8 ) ) then
         N_Dim_Htot=N_Dim_Htot+1    ! only the exchange is necessary
         cycle
      else
         N_Dim_Htot=N_Dim_Htot+indexNN(i,1)    ! more tensors are necessary to take into account the DMI symmetries
         cycle
      endif
   endif
   N_Dim_Htot=N_Dim_Htot+1

enddo

allocate(Hamiltonian%total_shell(N_max_shell+1))
do i=1,N_max_shell+1
   allocate(Hamiltonian%total_shell(i)%atom(N_max_nei))
   do j=1,N_max_nei
      allocate(Hamiltonian%total_shell(i)%atom(j)%H(3,3))
      Hamiltonian%total_shell(i)%atom(j)%H=0.0d0
   enddo
enddo

! set up the total Hamiltonian per site
! the Zeeman is not inside because it is a simple scalar produc in the total energy and in the Beff
if (N_coeff_ani.ne.0) Hamiltonian%total_shell(1)%atom(1)%H=Hamiltonian%total_shell(1)%atom(1)%H+Hamiltonian%ani(:,:,1)
do i=2,N_max_shell+1

   if (i-1.le.N_coeff_DMI) then

      if ( .not.all( abs(Hamiltonian%DMI(:,:,i-1)).lt.1.0d-8 ) ) then
         do l=1,indexNN(i-1,1)

            if (i-1.le.N_coeff_Exch) Hamiltonian%total_shell(i)%atom(l)%H=Hamiltonian%total_shell(i)%atom(l)%H+Hamiltonian%exchange(:,:,i-1)

            Hamiltonian%total_shell(i)%atom(l)%H(1,2)=Hamiltonian%total_shell(i)%atom(l)%H(1,2)+Hamiltonian%DMI(1,2,i-1)*D(l,3,1)
            Hamiltonian%total_shell(i)%atom(l)%H(1,3)=Hamiltonian%total_shell(i)%atom(l)%H(1,3)+Hamiltonian%DMI(1,3,i-1)*D(l,2,1)
            Hamiltonian%total_shell(i)%atom(l)%H(2,3)=Hamiltonian%total_shell(i)%atom(l)%H(2,3)+Hamiltonian%DMI(2,3,i-1)*D(l,1,1)

            Hamiltonian%total_shell(i)%atom(l)%H(2,1)=Hamiltonian%total_shell(i)%atom(l)%H(2,1)+Hamiltonian%DMI(2,1,i-1)*D(l,3,1)
            Hamiltonian%total_shell(i)%atom(l)%H(3,1)=Hamiltonian%total_shell(i)%atom(l)%H(3,1)+Hamiltonian%DMI(3,1,i-1)*D(l,2,1)
            Hamiltonian%total_shell(i)%atom(l)%H(3,2)=Hamiltonian%total_shell(i)%atom(l)%H(3,2)+Hamiltonian%DMI(3,2,i-1)*D(l,1,1)

         enddo
      endif
   else

      if (i-1.le.N_coeff_Exch) Hamiltonian%total_shell(i)%atom(1)%H=Hamiltonian%total_shell(i)%atom(1)%H+Hamiltonian%exchange(:,:,i-1)

   endif

   if (i.le.N_coeff_ani) Hamiltonian%total_shell(i)%atom(1)%H=Hamiltonian%total_shell(i)%atom(1)%H+Hamiltonian%ani(:,:,i)
enddo

#ifdef CPP_DEBUG
do i=1,N_max_shell+1
   write(*,*) 'shell', i

   do j=1,N_max_nei

      write(6,*) Hamiltonian%total_shell(i)%atom(j)%H(1,:)
      write(6,*) Hamiltonian%total_shell(i)%atom(j)%H(2,:)
      write(6,*) Hamiltonian%total_shell(i)%atom(j)%H(3,:)
   enddo
enddo
#endif

end subroutine get_total_Hamiltonian

!
! function that determines how many coefficients are none 0 in the exchange or the DMI or whatever
!
integer function number_nonzero_coeff_2d(coeff,name)
implicit none
real(kind=8), intent(in) :: coeff(:,:)
character(len=*), intent(in) :: name
! internal
integer :: i,N(2)
real(kind=8) :: norm

N=shape(coeff(:,:))
number_nonzero_coeff_2d=0

! loop over the number of coefficients
do i=1,N(2)
   norm=sqrt(coeff(1,i)**2+coeff(2,i)**2+coeff(3,i)**2)
   if (norm.gt.1.0d-8) number_nonzero_coeff_2d=i
enddo

if (number_nonzero_coeff_2d.eq.0) then
   write(6,'(2a)') 'no non-zero coefficients found for ',name
endif

end function number_nonzero_coeff_2d

integer function number_nonzero_coeff_1d(coeff,name)
implicit none
real(kind=8), intent(in) :: coeff(:)
character(len=*), intent(in) :: name
! internal
integer :: i,N

N=size(coeff)
number_nonzero_coeff_1d=0

! loop over the number of coefficients
do i=1,N
   if (abs(coeff(i)).gt.1.0d-8) number_nonzero_coeff_1d=i
enddo

if (number_nonzero_coeff_1d.eq.0) then
   write(6,'(2a)') 'no non-zero coefficients found for ',name
endif

end function number_nonzero_coeff_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! count the number of shell needed for the simulation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer function get_number_shell(Hamiltonian)
implicit none
type(Coeff_Ham), intent(in) :: Hamiltonian
! internal
integer :: N_coeff_ani,N_coeff_DMI,N_coeff_Exch,N_coeff

get_number_shell=0

! obtain the size of all matrices
N_coeff_Exch=size(Hamiltonian%exchange,3)
N_coeff_ani=size(Hamiltonian%ani,3)
N_coeff_DMI=size(Hamiltonian%DMI,3)

N_coeff=max(N_coeff_Exch,N_coeff_ani,N_coeff_DMI)

if (N_coeff.eq.0) then
   write(6,'(a)') 'WARNING!!!! No shells were found in the Hamiltonian (H=0)'
else
   write(6,'(I6,a/)') N_coeff,' shells were found in the Hamiltonian'
endif

get_number_shell=N_coeff

end function get_number_shell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! prepare the Hamiltonian matrix and the pointers for the matrix Hamiltonian static association
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine associate_energy_Hamiltonian(my_lattice,tableNN,indexNN)
use m_derived_types, only : Coeff_Ham,lattice
implicit none
type(lattice), intent(in) :: my_lattice
integer, intent(in) :: tableNN(:,:,:,:,:,:) !!tableNN(4,N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),count(my_motif%i_mom)
integer, intent(in) :: indexNN(:)
! internal
!type(operator_real),save :: energy
! slope of the sums
integer :: Nspin,all_size(4),shape_tableNN(6)

all_size=shape(my_lattice%l_modes)
Nspin=product(all_size)
shape_tableNN=shape(tableNN)

allocate(energy%value(Nspin,Nspin))
energy%nline=Nspin
energy%ncolumn=Nspin

call dissociate(energy,Nspin,Nspin)

call associate_pointer(energy,Hamiltonian%total_shell,my_lattice,tableNN,indexNN)

end subroutine associate_energy_Hamiltonian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! optimize the calculation time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine associate_line_Hamiltonian(N_cell,E_line,spin,mode_E_column)
implicit none
integer, intent(in) :: N_cell
type(vec_point),intent(in) :: spin(:)
type(point_shell_Operator), intent(inout) :: E_line(:)
type(point_shell_mode), intent(inout) :: mode_E_column(:)
!internal energy
integer :: i,k,j,n_atom_shell

do i=1,N_cell
   ! count the number of shell per line
   n_atom_shell=0
   do j=1,N_cell
      if (associated(energy%value(j,i)%Op_loc)) n_atom_shell=n_atom_shell+1
   enddo

   allocate(E_line(i)%shell(n_atom_shell))
   allocate(mode_E_column(i)%shell(n_atom_shell))

   ! the k=1 case should be the case i=j
   ! it should be the first atom in the shell
   k=1
   do j=1,N_cell

      if (associated(energy%value(j,i)%Op_loc)) then
         ! taking care of the on-site term
         if (i.eq.j) then
            E_line(i)%shell(1)%Op_loc=>energy%value(j,i)%Op_loc
            mode_E_column(i)%shell(1)%w=>spin(j)%w
            cycle
         endif

         k=k+1
         E_line(i)%shell(k)%Op_loc=>energy%value(j,i)%Op_loc
         mode_E_column(i)%shell(k)%w=>spin(j)%w
      endif

   enddo

   if (k.ne.n_atom_shell) stop 'error in associate_line_Hamiltonian'

enddo

end subroutine associate_line_Hamiltonian

end module m_energy_commons
