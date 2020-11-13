module m_setup_simu
implicit none

contains
subroutine setup_simu(io_simu,my_lattice,my_motif,ext_param,Ham)
use m_derived_types
use m_energyfield, only : init_Energy_distrib
use m_energy_commons
use m_internal_fields_commons
use m_fftw
use m_lattice
use m_setup_DM
use m_user_info
use m_sym_utils
use m_table_dist
use m_mapping
use m_indexation
use m_arrange_neigh
use m_topoplot
use m_get_position
use m_external_fields
use m_excitations
use m_io_files_utils
use m_io_utils
use m_dipolar_field
use m_null
use m_rw_TB, only : rw_TB, check_activate_TB, get_nb_orbitals
use m_summer_exp, only : get_coeff_TStra
use m_input_H_types

use m_rw_H
use m_H_public

#ifdef CPP_MPI
      use m_make_box
      use m_split_work
      use m_mpi_prop, only : isize,irank,irank_working,N,start
#endif

! this subroutine is used only to setup the simulation box
! it reads first the parameters of the simulation i.e. inp file
! then it reads the lattice
! this order aims at not taking care of too many neigbours if the Jij are set to 0
type(io_parameter), intent(out) :: io_simu
type(lattice), intent(inout) :: my_lattice
type(t_cell), intent(out) :: my_motif
type(simulation_parameters),intent (inout) :: ext_param
class(t_H),intent(inout),allocatable      ::  Ham(:)
! variable of the system
real(kind=8), allocatable :: tabledist(:,:),DM_vector(:,:,:)
integer, allocatable :: indexNN(:,:),tableNN(:,:,:,:,:,:)
real (kind=8), allocatable :: pos(:,:,:,:,:)
integer :: tot_N_Nneigh,io
real(kind=8) :: time
! dummy variable
integer :: dim_lat(3),n_mag,n_DMI,N_Nneigh
!checking various files
logical :: i_usestruct
! check the allocation of memory
integer :: alloc_check
! Nb orbitals in TB
integer :: nb_orbitals = 0
! Hamiltonian input parameters
type(io_h)  ::  H_io


! innitialisation of dummy
i_usestruct=.False.
alloc_check=0
time=0.0d0

call user_info(6,time,'entering the setup routine',.True.)

! read the important inputs

time=0.0d0
call user_info(6,time,'reading the lattice in the input file',.False.)

call rw_lattice(my_lattice)

! read the important inputs
call user_info(6,time,'reading the motif in the input file',.False.)

call rw_motif(my_motif,my_lattice)

!read the Hamiltonian parameters
Call rw_H(H_io)


! read the TB parameters
call rw_TB('input')
if (check_activate_TB()) nb_orbitals=get_nb_orbitals()

! find the symmetry of the lattice here

call user_info(6,time,'reading the input file',.false.)
time=0.0d0

call inp_rw(io_simu)

call user_info(6,time,'Read external parameters',.false.)
time=0.0d0

call ext_param_rw(ext_param)

call user_info(6,time,'reading the Hamiltonian in the input file',.false.)

!call get_excitations('input')

call user_info(6,time,'done',.true.)

#ifdef CPP_MPI
  if (irank.eq.0) write(6,'(a)') 'checking for the presence of electric field'
#else
  write(6,'(a,/)') 'checking for the presence of electric field'
#endif

! check for electric field
call rw_efield(my_lattice%dim_lat,my_lattice%areal)
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! create the lattice depending on the simulation that one chooses
call user_info(6,time,'allocating the spin, table of neighbors and the index...',.false.)
time=0.0d0

call create_lattice(my_lattice,my_motif,ext_param,nb_orbitals)

!stop 'toto'
! note for Bertrand: Continue from here
dim_lat=my_lattice%dim_lat
n_mag=count(my_motif%atomic(:)%moment.gt.0.0d0)

! check if you put the Strasbourg temperature
call get_coeff_TStra(my_lattice,my_motif)

! setup the Hamiltonian of the system
! One need the dimension of the order parameter before the hamiltonians can be created
call get_Hamiltonians('input',my_motif%atomic(1)%moment,my_lattice%dim_mode)


! get the table of neighbors and the numbers of atom per shell
N_Nneigh=H_io%get_shell_number()
allocate(tabledist(N_Nneigh,1),stat=alloc_check)
if (alloc_check.ne.0) write(6,'(a)') 'out of memory cannot allocate table of distance'
tabledist=0.0d0
allocate(indexNN(N_Nneigh,1),stat=alloc_check)
if (alloc_check.ne.0) write(6,'(a)') 'out of memory cannot allocate indexNN'
indexNN=0
call get_table_of_distance(my_lattice%areal,N_Nneigh,my_lattice%world,my_motif,n_mag,tabledist)
call get_num_neighbors(N_Nneigh,tabledist,my_lattice%areal,my_lattice%world,my_motif,indexNN)

! prepare the lattice
call user_info(6,time,'initializing the spin structure',.false.)

call InitSpin(my_lattice,my_motif,ext_param)

call user_info(6,time,'done',.false.)

! allocate table of neighbors and masque
tot_N_Nneigh=sum(indexNN(1:N_Nneigh,1),1)

allocate(tableNN(5,tot_N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),n_mag),stat=alloc_check)
if (alloc_check.ne.0) write(6,'(a)') 'out of memory cannot allocate tableNN'

tableNN=0

!-------------------------------------------------
!mapping of the neighbours
call user_info(6,time,'Calculating the table of neighbors (this can take time)',.true.)

allocate(pos(3,dim_lat(1),dim_lat(2),dim_lat(3),n_mag),stat=alloc_check)
if (alloc_check.ne.0) write(6,'(a)') 'out of memory cannot allocate position for the mapping procedure'

! get position
pos=0.0d0
call get_position(pos,dim_lat,my_lattice%areal,my_motif)
! write the positions into a file
io=open_file_write('positions.dat')
call dump_config(io,pos)
call close_file('positions.dat',io)

! get table of neighbors
call mapping(tabledist,N_Nneigh,my_motif,indexNN,tableNN,my_lattice)

call user_info(6,time,'done',.true.)

! prepare the external electromagnetic fields
call user_info(6,time,'get the external magnetic field matrix',.false.)

call get_EM_external_fields(ext_param%H_ext%value,ext_param%E_ext%value,my_motif)

call initialize_external_fields(my_lattice)

call user_info(6,time,'done',.false.)

! setup the different parameters for the interaction DM, 4-spin... 
! hard part: setup DM interaction
! I send an example neighbor table to setup the DM

call user_info(6,time,'Calculating the DM vectors for each shells',.false.)

n_DMI=get_number_DMI('input')

if (n_DMI.ne.0) then
   write(6,'(I4,a)') n_DMI,' DMI found'

   write(*,*) 'number of first nearest neighbor', sum(indexNN(1:n_DMI,1))
   allocate(DM_vector(sum(indexNN(1:n_DMI,1)),3,1))

   DM_vector=0.0d0
   call setup_DM_vector(indexNN,n_DMI,my_lattice,my_motif,DM_vector,tabledist)

   call user_info(6,time,'done',.true.)

! the DM and the neighbors have to turn in the same direction. therefore the matrix of the neighbors
! have to be rearranged
    call user_info(6,time,'Re-aranging the position of the DM vectors',.false.)

    call arrange_DM(DM_vector,n_DMI)
!    call arrange_neigh(DM_vector,tableNN,indexNN,my_lattice%dim_lat,my_lattice%areal,n_DMI)

    call user_info(6,time,'done',.true.)
else
    allocate(DM_vector(1,1,1))
    DM_vector=0.0d0
endif

!if (Hamiltonian%i_four) then
!     call user_info(6,time,'Calculating the position of the corners of the unit cell for the 4-spin',.false.)
!
!       call givemecorners(my_lattice%areal,abs(order_zaxis(my_lattice%areal)))
!
!     call user_info(6,time,'done',.true.)
!else
       allocate(corners(1,1))
!endif

! do a first FFT for the initial magnetic configuration
!call fft(my_lattice,my_motif)



!
! get the null matrix in case no boundary conditions
!
call get_null_matrix(my_lattice%dim_mode)

!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!
!!!! important part that associates the different NxN matrices with the local static operators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!call user_info(6,time,'associating the operator matrices with the different local variables',.true.)
!
! prepare the total Hamiltonian depending on the different interactions
!
! call  get_total_Hamiltonian(DM_vector,indexNN)
!
! update the Hamiltonian to take into account the lattice i.e. SOC (DMI...)
!
!
!call associate_energy_Hamiltonian(my_lattice,tableNN,indexNN(:,1))
!
!call associate_internal_Beff(my_lattice,tableNN,indexNN(:,1))
!
!call associate_external_fields(tableNN)
!
!call user_info(6,time,'done',.true.)

!!--------------------------------------------------------------------

!call user_info(6,time,'allocating the different variables for the energy/angle... density distribution',.true.)
!
!if (io_simu%io_Energy_Distrib) call init_Energy_distrib(my_lattice,tableNN,indexNN(:,1))
!
!call user_info(6,time,'done',.true.)





!
! Check the presence of the dipole dipole
!
!
call get_ham_dipole('input',my_lattice,my_motif)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! CREATE NEW HAMILTONIAN TYPE
!old
!Call get_Htype(Ham)
!Call Ham%set_H(energy,my_lattice)
!
Call set_Hamiltonians(Ham,H_io,tableNN,indexNN(:,1),DM_vector,my_lattice)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
deallocate(tabledist,tableNN,indexNN)


#ifdef CPP_MPI
if (irank.eq.0) write(6,'(/,a/)') 'the setup of the simulation is over'
#else
write(6,'(/,a,/)') 'the setup of the simulation is over'
#endif

if (io_simu%io_fft_Xstruct) call set_k_mesh('input',my_lattice)


!!!!!!!!!!!!!! end of the setup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine setup_simu

subroutine set_Hamiltonians(Ham,H_io,tableNN,indexNN,DM_vector,lat)
    use m_derived_types
    use m_H_public
    use m_input_H_types 
    use m_get_position,only :get_position_ND_to_1D 
    
    use m_symmetry_operators
    use m_anisotropy_heisenberg,only: get_anisotropy_H
    !use m_exchange_heisenberg,only: get_exchange_H,exchange
    use m_zeeman,only: get_zeeman_H
    use m_exchange_heisenberg_J, only: get_exchange_J
    use m_exchange_heisenberg_D, only: get_exchange_D
    use m_coupling_ME_J,only: get_coupling_ME_J
    use m_coupling_ME_D,only: get_coupling_ME_D
    class(t_H),allocatable,intent(out)  :: Ham(:)
    type(io_h),intent(in)               :: H_io
    real(8), intent(in) :: DM_vector(:,:,:)
    integer, intent(in) :: tableNN(:,:,:,:,:,:) !!tableNN(4,N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),count(my_motif%i_mom)
    integer, intent(in) :: indexNN(:)
    type(lattice), intent(inout) :: lat

    integer :: i_H,N_ham
    logical :: use_Ham(6)

    use_ham(1)=H_io%J%is_set
    use_ham(2)=H_io%D%is_set
    use_ham(3)=H_io%aniso%is_set
    use_ham(4)=H_io%zeeman%is_set
    use_ham(5)=H_io%ME_J%is_set
    use_ham(6)=H_io%ME_D%is_set
    N_ham=count(use_ham)
    Call get_Htype_N(Ham,N_ham)

    i_H=1 
    !exchange_J (without DMI)
    if(use_ham(1))then
        Call get_exchange_J(Ham(i_H),H_io%J,tableNN,indexNN,lat)
        if(Ham(i_H)%is_set()) i_H=i_H+1
    endif
    !exchange_D (only DMI)
    if(use_ham(2))then
        Call get_exchange_D(Ham(i_H),H_io%D,tableNN,indexNN,lat,DM_vector)
        if(Ham(i_H)%is_set()) i_H=i_H+1
    endif
    !anisotropy
    if(use_ham(3))then
        Call get_anisotropy_H(Ham(i_H),H_io%aniso,lat)
        if(Ham(i_H)%is_set()) i_H=i_H+1
    endif
    !zeeman
    if(use_ham(4))then
        Call get_zeeman_H(Ham(i_H),H_io%zeeman,lat)
        if(Ham(i_H)%is_set()) i_H=i_H+1
    endif
    !ME-coupling symmetric (J)
    if(use_ham(5))then
        Call get_coupling_ME_J(Ham(i_H),H_io%ME_J,tableNN,indexNN,lat)
        if(Ham(i_H)%is_set()) i_H=i_H+1
    endif
    !ME-coupling antisymmetric (D)
    if(use_ham(5))then
        Call get_coupling_ME_D(Ham(i_H),H_io%ME_D,tableNN,indexNN,lat)
        if(Ham(i_H)%is_set()) i_H=i_H+1
    endif

    do i_H=1,N_ham
        if(.not. Ham(i_h)%is_set()) STOP "not all Hamiltonians are set"
        Call Ham(i_h)%optimize()
    enddo


end subroutine


end module
