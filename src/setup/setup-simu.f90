module m_setup_simu
implicit none

contains
subroutine setup_simu(io_simu,my_lattice,ext_param,Ham_res,Ham_comb,H_res,H_comb)
    !main setup routine (right now should only be called from one MPI-thread, bcast afterwards)
    use m_derived_types, only: io_parameter,t_cell,simulation_parameters
    use m_fftw, only: set_k_mesh
    use m_setup_DM
    use m_user_info
    use m_table_dist
    use m_mapping
    use m_indexation
    use m_arrange_neigh
    use m_get_position
    use m_io_files_utils, only: open_file_write, close_file
    use m_io_utils, only: dump_config
    use m_set_Hamiltonians,only: set_Hamiltonians
    use m_set_fft_Hamiltonians,only: set_fft_Hamiltonians
    use m_rw_extpar, only: extpar_input, rw_extpar
    use m_orders_initialize, only: orders_initialize 
    use m_rw_cell
    use m_rw_H
    use m_H_public
    use m_fft_H_public
    use m_neighbor_type
    use m_hamiltonian_collection, only: hamiltonian
    use m_lattice_position, only: get_pos_mag, print_pos_mag, print_pos_atom
    
    ! this subroutine is used only to setup the simulation box
    ! it reads first the parameters of the simulation i.e. inp file
    ! then it reads the lattice
    ! this order aims at not taking care of too many neigbours if the Jij are set to 0
    type(io_parameter), intent(out)             :: io_simu
    type(lattice), intent(inout)                :: my_lattice
    type(simulation_parameters),intent (inout)  :: ext_param
    class(t_H),intent(inout),allocatable        :: Ham_res(:), Ham_comb(:)
    type(hamiltonian),intent(inout)             :: H_res,H_comb

    class(fft_H),allocatable    :: ffT_Ham_res(:)
    class(fft_H),allocatable    :: fft_Ham_comb(:)

    ! variable of the system
    type(t_cell)        :: my_motif
    real(8),allocatable :: pos(:)
    integer             :: io
    real(8)             :: time
    type(extpar_input)  :: extpar_io
    ! dummy variable
    logical             :: keep_resolved_H
    ! check the allocation of memory
    ! Hamiltonian input parameters
    type(io_h)          :: H_io
    
    !initialisation of dummy
    time=0.0d0
    call user_info(6,time,'entering the setup routine',.True.)
    
    ! read the important inputs
    time=0.0d0
    call user_info(6,time,'reading the lattice in the input file',.False.)
    call rw_lattice(my_lattice)
    
    ! read the important inputs
    call user_info(6,time,'reading the motif in the input file',.False.)
    Call read_cell(my_motif,my_lattice)

    !read the Hamiltonian parameters
    Call rw_H(H_io)
    
    ! find the symmetry of the lattice here
    call user_info(6,time,'reading the input file',.false.)
    time=0.0d0
    call inp_rw(io_simu)
    
    call user_info(6,time,'Read external parameters',.false.)
    time=0.0d0
    call ext_param_rw(ext_param)
    Call rw_extpar(extpar_io)
    write(6,'(a,/)') 'checking for the presence of electric field'
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! create the lattice depending on the simulation that one chooses
    call user_info(6,time,'allocating the spin, table of neighbors and the index...',.false.)
    time=0.0d0
    
    Call my_lattice%init_order(my_motif,extpar_io)
    
    !-------------------------------------------------
    
    ! prepare the lattice
    write(*,'(/)')
    call user_info(6,time,'Start initializing the order parameters',.false.)
    Call orders_initialize(my_lattice,extpar_io)
    call user_info(6,time,'Finished initializing the order parameters',.false.)
    write(*,'(/)')

    ! write the positions into a file
    Call print_pos_mag (my_lattice,'positions.dat')
    Call print_pos_atom(my_lattice,'position_atoms.dat')


    write(6,'(///)') 
    call user_info(6,time,'Start setting Hamiltonians',.false.)
    keep_resolved_H=io_simu%io_Energy_detail.or..True.  !need to change where all Hamiltonian data is kept
    Call set_Ham_mode_io()
    Call set_Hamiltonians(Ham_res,Ham_comb,keep_resolved_H,H_io,my_lattice)
    call user_info(6,time,'finished setting Hamiltonians',.false.)

    Call set_fft_H_mode_io()
    call user_info(6,time,'Start setting fft-Hamiltonians',.false.)
    Call set_fft_Hamiltonians(fft_Ham_res,fft_Ham_comb,keep_resolved_H,H_io,my_lattice)
    call user_info(6,time,'finished setting fft-Hamiltonians',.false.)


    Call H_comb%init_H_cp(my_lattice,Ham_comb,fft_Ham_comb)   !later change to move if the hamiltonian-type array is no longer necessary in main routine
    if(keep_resolved_H) Call H_res%init_H_cp(my_lattice,Ham_res,fft_Ham_res)   !later change to move if the hamiltonian-type array is no longer necessary in main routine

    if (io_simu%io_fft_Xstruct) call set_k_mesh('input',my_lattice)
    
    write(6,'(/,a,/)') 'the setup of the simulation is over'
    write(6,'(I6,a)') my_lattice%ncell, ' unit cells'
    write(6,'(a)') '-----------------------------------------------'
    write(6,'(a)') ''
    write(6,'(a)') '-----------------------------------------------'
end subroutine setup_simu

end module
