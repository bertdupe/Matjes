subroutine setup_simu(my_simu,io_simu,my_lattice,my_motif,ext_param,state)
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
use mtprng
use m_total_energy
use m_get_position
use m_external_fields
use m_excitations
use m_io_files_utils
use m_io_utils
use m_dipolar_field
use m_null

#ifdef CPP_MPI
      use m_make_box
      use m_split_work
      use m_mpi_prop, only : isize,irank,irank_working,N,start
#endif

implicit none
! this subroutine is used only to setup the simulation box
! it reads first the parameters of the simulation i.e. inp file
! then it reads the lattice
! this order aims at not taking care of too many neigbours if the Jij are set to 0
type(io_parameter), intent(out) :: io_simu
type(bool_var), intent(in) :: my_simu
type(lattice), intent(inout) :: my_lattice
type(cell), intent(out) :: my_motif
type (mtprng_state),intent(inout) :: state
type(simulation_parameters),intent (inout) :: ext_param
! variable of the system
real(kind=8), allocatable :: tabledist(:,:),DM_vector(:,:,:)
integer, allocatable :: indexNN(:,:),tableNN(:,:,:,:,:,:),masque(:,:,:,:)
real (kind=8), allocatable :: pos(:,:,:,:,:)
integer :: tot_N_Nneigh,io,dim_ham
real(kind=8) :: time
! dummy variable
integer :: i,dim_lat(3),n_mag,n_DMI,N_Nneigh
!checking various files
logical :: topoonly,i_exi,i_usestruct
! check the allocation of memory
integer :: alloc_check

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

call create_lattice(my_lattice,my_motif,ext_param)

dim_lat=my_lattice%dim_lat
n_mag=count(my_motif%atomic(:)%moment.gt.0.0d0)

! setup the Hamiltonian of the system
! One need the dimension of the order parameter before the hamiltonians can be created
call get_Hamiltonians('input',my_motif%atomic(1)%moment,my_lattice%dim_mode)

N_Nneigh=get_number_shell()

allocate(indexNN(N_Nneigh,1),stat=alloc_check)
if (alloc_check.ne.0) write(6,'(a)') 'out of memory cannot allocate indexNN'
allocate(tabledist(N_Nneigh,1),stat=alloc_check)
if (alloc_check.ne.0) write(6,'(a)') 'out of memory cannot allocate table of distance'
indexNN=0
tabledist=0.0d0

call get_table_of_distance(my_lattice%areal,N_Nneigh,my_lattice%world,my_motif,indexNN,n_mag,tabledist)

call get_num_neighbors(N_Nneigh,tabledist,my_lattice%areal,my_lattice%world,my_motif,indexNN)

! prepare the lattice
call user_info(6,time,'initializing the spin structure',.false.)

call InitSpin(my_lattice,my_motif,state,ext_param)

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
call mapping(tabledist,N_Nneigh,my_motif,indexNN,tableNN,my_lattice,pos)

call user_info(6,time,'done',.true.)

! prepare the external electromagnetic fields
call user_info(6,time,'get the external magnetic field matrix',.false.)

call get_EM_external_fields(ext_param%H_ext%value,ext_param%E_ext%value,my_motif)

call initialize_external_fields(my_lattice)

call user_info(6,time,'done',.false.)

!check for the structure.xyz file for shape modification
inquire (file='structure.xyz',exist=i_usestruct)
if (i_usestruct) call user_def_struct(masque,pos, &
     & my_lattice%dim_lat,count(my_motif%i_mom),tot_N_Nneigh+1,my_lattice%areal)

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
   call setup_DM_vector(indexNN,n_DMI,my_lattice,my_motif,DM_vector)

   call user_info(6,time,'done',.true.)

! the DM and the neighbors have to turn in the same direction. therefore the matrix of the neighbors
! have to be rearranged
    call user_info(6,time,'Re-aranging the position of the DM vectors',.false.)

    if (size(my_lattice%world).eq.1) then

    elseif (size(my_lattice%world).eq.2) then

         call arrange_neigh(DM_vector(:,:,1),tableNN(:,:,:,:,1,:),indexNN(:,1),my_lattice%dim_lat,my_lattice%areal)

!      else
!         call arrange_neigh(DM_vector,tableNN(:,:,:,:,1,:),indexNN(:,:),my_lattice%dim_lat,my_lattice%areal,masque)

    else
      write(*,*) 'not coded'
    endif

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

! setup the energy operator
call user_info(6,time,'dealing with the z-direction',.false.)

!! check the z direction structure
!      call setup_zdir(phase,tot_N_Nneigh,motif)

call user_info(6,time,'done',.true.)

! do a first FFT for the initial magnetic configuration
!call fft(my_lattice,my_motif)



!
! get the null matrix in case no boundary conditions
!
call get_null_matrix(my_lattice%dim_mode)

!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!
!!!! important part that associates the different NxN matrices with the local static operators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call user_info(6,time,'associating the operator matrices with the different local variables',.true.)
!
! prepare the total Hamiltonian depending on the different interactions
!
call  get_total_Hamiltonian(DM_vector,indexNN)
!
! update the Hamiltonian to take into account the lattice i.e. SOC (DMI...)
!

call associate_energy_Hamiltonian(my_lattice,tableNN,indexNN(:,1))

call associate_internal_Beff(my_lattice,tableNN,indexNN(:,1))

call associate_external_fields(tableNN)

call user_info(6,time,'done',.true.)

!!--------------------------------------------------------------------

call user_info(6,time,'allocating the different variables for the energy/angle... density distribution',.true.)

if (io_simu%io_Energy_Distrib) call init_Energy_distrib(my_lattice,tableNN,indexNN(:,1))

call user_info(6,time,'done',.true.)

!c Calculation of the dispertion in the full BZ
! this is here because we have the table of distances
!      if (dispersion) then
!       call fullBZ(tabledist,N_Nneigh,Nei_z,phase)
!      endif

!deallocate(tabledist,tableNN,indexNN)



!
! Check the presence of the dipole dipole
!
!
call get_ham_dipole('input',my_lattice,my_motif)

#ifdef CPP_MPI
if (irank.eq.0) write(6,'(/,a/)') 'the setup of the simulation is over'
#else
write(6,'(/,a,/)') 'the setup of the simulation is over'
#endif

if (io_simu%io_fft_Xstruct) call get_k_mesh('input',my_lattice)

!!!!!!!!!!!!!! end of the setup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine setup_simu
