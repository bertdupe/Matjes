! format des différentes variables
!type lattice
!     real(kind=8) :: areal(3,3),astar(3,3),alat(3)
!     integer :: dim_lat(3),n_system,dim_mode
!     integer, allocatable :: world(:)
!     logical :: boundary(3)
!! Table of pointer
!     type(vec_point), allocatable :: l_modes(:,:,:,:)
!end type lattice





subroutine tightbinding(my_lattice,my_motif,io_simu,ext_param)
    use m_basic_types, only : vec_point
    use m_derived_types, only : cell,lattice,io_parameter,simulation_parameters
    use m_local_energy, only : get_E_matrix
    use m_fftw
    use m_lattice
    use m_operator_pointer_utils
    use m_kmesh

    implicit none
    ! internal parameter
    type(io_parameter), intent(in) :: io_simu
    type(lattice), intent(in) :: my_lattice
    type(cell), intent(in) :: my_motif
    type(simulation_parameters), intent(in) :: ext_param

    ! N_cell is the variable that will contain the number of unit
    ! cells in the simulation
    integer :: N_cell
    ! shape_lattice is an array of 4 integer that will contain the
    ! diemnsions along x, y, and z, and the number of atoms per unit
    ! cell (for exemple, 1 in BCC, 4 in FCC, etc.)
    integer :: shape_lattice(4)
    ! all_mode is an array of custom type vec_point that will
    ! contain all the modes present in the simulation
    type(vec_point),allocatable :: all_mode(:)
    ! dim_mode contains the dimension of the mode array
    integer :: dimm

    call get_k_mesh('input',my_lattice)

    shape_lattice=shape(my_lattice%l_modes)
    N_cell=product(shape_lattice)

    allocate(all_mode(N_cell))

    call associate_pointer(all_mode,my_lattice)
    call get_E_matrix(my_lattice%dim_mode)

    ! Perform the FFT of the Hamiltonian
    ! In that function, the arguments are:
    !   _ kmesh: array of kpoints
    !   _ -1.0: the direction of the transform (-1.0=direct, +1.0=reverse)
    !   _ my_motif%pos: position of the atoms in the lattice
    !   _ field: ???
    !   _ Nsize: number of sites in the lattice (here equal to N_cell ?)
    !   _ my_lattice%dim_mode: dimension of the order parameter
    call get_FFT_vec_point(kmesh,-1.0,my_motif%pos,all_mode,N_cell,my_lattice%dim_mode)

    write(*,*) 'Coucou from file ', __FILE__, ' at line ', __LINE__
end subroutine tightbinding