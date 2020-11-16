module m_setup_simu
implicit none

contains
subroutine setup_simu(io_simu,my_lattice,my_motif,ext_param,Ham_res,Ham_comb)
    use m_derived_types
    !use m_energyfield, only : init_Energy_distrib
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
    use m_set_Hamiltonians,only: set_Hamiltonians
    
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
    class(t_H),intent(inout),allocatable      ::  Ham_res(:), Ham_comb(:)
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
    
    dim_lat=my_lattice%dim_lat
    n_mag=count(my_motif%atomic(:)%moment.gt.0.0d0)
    
    ! check if you put the Strasbourg temperature
    call get_coeff_TStra(my_lattice,my_motif)
    
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
    Call user_info(6,time,'get the external magnetic field matrix',.false.)
    Call get_EM_external_fields(ext_param%H_ext%value,ext_param%E_ext%value,my_motif)
    Call initialize_external_fields(my_lattice)
    Call user_info(6,time,'done',.false.)
    
    ! setup the different parameters for the interaction DM, 4-spin... 
    ! hard part: setup DM interaction
    ! I send an example neighbor table to setup the DM
    call user_info(6,time,'Calculating the DM vectors for each shells',.false.)
    n_DMI=count(H_io%D%val/=0.0d0)
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
        call user_info(6,time,'done',.true.)
    else
        allocate(DM_vector(1,1,1))
        DM_vector=0.0d0
    endif
    
    ! Check the presence of the dipole dipole
    !call get_ham_dipole('input',my_lattice,my_motif) !DIPOLE HAS THE BE REIMPLEMENTED
    
    Call set_Hamiltonians(Ham_res,Ham_comb,io_simu%io_Energy_detail,H_io,tableNN,indexNN(:,1),DM_vector,my_lattice)
    
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

end module
