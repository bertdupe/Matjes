module m_setup_simu
use m_type_lattice, only: lattice
use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit, output_unit
use m_user_info
implicit none
private
public :: setup_simu

contains
subroutine setup_simu(io_simu,my_lattice,ext_param,Ham_res,Ham_comb,H_res,H_comb)
    !main setup routine (right now should only be called from one MPI-thread, bcast afterwards)
    use m_derived_types, only: io_parameter,t_cell,simulation_parameters
    use m_fftw, only: set_k_mesh
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
    use m_diagonalization_Hk
    use m_grp_sym
    use m_precision, only : set_EPS
    
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
    real(8)             :: time
    type(extpar_input)  :: extpar_io
    ! dummy variable
    logical             :: keep_resolved_H
    ! check the allocation of memory
    ! Hamiltonian input parameters
    type(io_h)          :: H_io
    
    ! call for a truncation to set the 0
    call set_EPS()

    !initialisation of dummy
    time=0.0d0
    call user_info(output_unit,time,'entering the setup routine',.True.)
    
    ! read the important inputs
    time=0.0d0
    call user_info(output_unit,time,'reading the lattice in the input file',.False.)
    call rw_lattice(my_lattice)
    
    ! read the important inputs
    call user_info(output_unit,time,'reading the motif in the input file',.False.)
    Call read_cell(my_motif,my_lattice)

    !read the Hamiltonian parameters
    Call rw_H(H_io)
    
    ! find the symmetry of the lattice here
    call user_info(output_unit,time,'reading the input file',.false.)
    time=0.0d0
    call inp_rw(io_simu)
    
    call user_info(output_unit,time,'Read external parameters',.false.)
    time=0.0d0
    Call rw_extpar(extpar_io,my_lattice%areal)
    Call ext_param%set(extpar_io%H,extpar_io%E,extpar_io%T) !set legacy external parameters from newer input (extpar_io)
    write(output_unit,'(a,/)') 'checking for the presence of electric field'
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! create the lattice depending on the simulation that one chooses
    call user_info(output_unit,time,'allocating the spin, table of neighbors and the index...',.false.)
    time=0.0d0
    
    Call my_lattice%init_order(my_motif,extpar_io)
    
    !-------------------------------------------------
    
    ! prepare the lattice
    call user_info(output_unit,time,'Start initializing the order parameters',.false.)
    Call orders_initialize(my_lattice,extpar_io)
    call user_info(output_unit,time,'Finished initializing the order parameters',.false.)

    ! write the positions into a file
    Call print_positions(my_lattice,time)

    ! get the space group and the point group
    call find_group(my_lattice,my_motif)

    write(output_unit,'(///)') 
    call user_info(output_unit,time,'Start setting Hamiltonians',.false.)
    keep_resolved_H=io_simu%io_Energy_detail.or..True.  !need to change where all Hamiltonian data is kept
    Call set_Ham_mode_io()
    Call set_Hamiltonians(Ham_res,Ham_comb,keep_resolved_H,H_io,my_lattice)
    call user_info(output_unit,time,'finished setting Hamiltonians',.false.)

    Call set_fft_H_mode_io()
    call user_info(output_unit,time,'Start setting fft-Hamiltonians for real space dynamics',.false.)
    Call set_fft_Hamiltonians(fft_Ham_res,fft_Ham_comb,keep_resolved_H,H_io,my_lattice)
    call user_info(output_unit,time,'finished setting fft-Hamiltonians',.false.)


    Call H_comb%init_H_cp(my_lattice,Ham_comb,fft_Ham_comb)   !later change to move if the hamiltonian-type array is no longer necessary in main routine
    if(keep_resolved_H) Call H_res%init_H_cp(my_lattice,Ham_res,fft_Ham_res)   !later change to move if the hamiltonian-type array is no longer necessary in main routine

    if (io_simu%io_fft_Xstruct) call set_k_mesh('input',my_lattice)
    

    if (io_simu%io_fft_Ham) then
       call user_info(6,time,'Start setting fft-Hamiltonians for diagonalization',.false.)
       call diagonalize_Ham_FT(H_io,my_lattice)
       call user_info(6,time,'End diagonalization',.false.)
    endif

    write(6,'(/,a,/)') 'the setup of the simulation is over'
    write(6,'(I15,a)') my_lattice%ncell, ' unit cells'
    write(6,'(a)') '-----------------------------------------------'
    write(6,'(a)') ''
    write(6,'(a)') '-----------------------------------------------'
end subroutine setup_simu

subroutine print_positions(lat,time)
    !print positions of atoms to file
    use m_lattice_position, only:  print_pos_ind
    type(lattice), intent(in)   :: lat 
    real(8),intent(inout)       :: time
    integer             :: i
    integer,allocatable :: ind(:)   !save atom indices (in unit-cell) considered in respective file
    integer             :: Nat      !number of atoms

    call user_info(output_unit,time,'Start printing atomic positions.',.false.)

    !print position of all atoms
    Nat=size(lat%cell%atomic)
    allocate(ind,source=[(i,i=1,Nat)])
    Call print_pos_ind(lat,ind,'position_atoms.dat')   !print position of all atoms to file
    deallocate(ind)

    !print position of all magnetic atoms
    if(lat%nmag>0)then
        Call lat%cell%ind_mag_all(ind)
        Call print_pos_ind (lat, ind, 'positions_magnetic.dat')
        deallocate(ind)
    endif

    !print position of all phonons
    if(lat%nph>0)then
        Call lat%cell%ind_Z_all(ind)
        Call print_pos_ind (lat, ind, 'positions_phonons.dat')
        deallocate(ind)
    endif
    call user_info(output_unit,time,'Finished printing atomic positions',.false.)
end subroutine
end module
