! format des diffÃ©rentes variables
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
    use  m_get_position

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
    ! dimensions along x, y, and z, and the number of atoms per unit
    ! cell (for exemple, 1 in BCC, 4 in FCC, etc.)
    integer :: shape_lattice(4)
    ! all_mode is an array of custom type vec_point that will
    ! contain all the modes present in the simulation
    type(vec_point),allocatable :: all_mode(:)
    ! sense gives the sense of the FFT transform
    real(kind=8) :: sense
    ! array containing all the positions in the lattice
    real(kind=8), allocatable :: all_positions(:,:,:,:,:)
    ! array that will contain the FFT coefficients
    complex(kind=16), allocatable :: array_FFT(:,:)
    complex(kind=16) :: testCOMPLEX(17,17)

    shape_lattice=shape(my_lattice%l_modes)
    N_cell=product(shape_lattice)

    allocate(all_mode(N_cell))

    ! We have access to the variable 'N_kpoint' because it is in the module
    ! m_fftw
    allocate(array_FFT(product(N_kpoint),product(N_kpoint)))

    allocate(all_positions(3,my_lattice%dim_lat(1),my_lattice%dim_lat(2),my_lattice%dim_lat(3),size(my_motif%atomic)))
    call get_position(all_positions,my_lattice%dim_lat,my_lattice%areal,my_motif)

    call associate_pointer(all_mode,my_lattice)
    call get_E_matrix(my_lattice%dim_mode)

    call get_k_mesh('input',my_lattice)

    sense=-1.0

    !call calcultae_fft(all_mode,all_positions,sense,array_FFT)
    !                                                         1
    !Error: There is no specific subroutine for the generic calculate_fft at (1)
    !call calculate_fft(all_mode,all_positions,sense,array_FFT)

    !undefined reference to calculate_fft_matrix
    call calculate_fft_matrix(all_mode,all_positions,sense,array_FFT)

    write(*,*) 'Coucou from file ', __FILE__, ' at line ', __LINE__
    write(*,*) 'shape_lattice = ', shape_lattice
end subroutine tightbinding