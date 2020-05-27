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
    use m_bandstructure
    use m_operator_pointer_utils
    use  m_get_position
    use m_io_utils
    use m_io_files_utils
    use m_DOS
    use m_wavefunction
    use m_energy_k

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
    real(kind=8), allocatable :: start_positions(:,:,:,:,:),pos(:,:),distances(:,:)

    complex(kind=16), allocatable :: dispersion(:), input_energy(:), DOS(:)
    integer :: io, i, nb_kpoints

    shape_lattice=shape(my_lattice%l_modes)
    N_cell=product(shape_lattice)

    call get_k_mesh('input',my_lattice)
    allocate( all_mode(N_cell), pos(3,N_cell), distances(3,N_cell) )
    pos=0.0d0
    distances=0.0d0

    ! We have access to the variable 'N_kpoint' because it is in the module
    ! m_fftw
    nb_kpoints = product(N_kpoint)
    allocate( dispersion(nb_kpoints) )
    dispersion=0.0d0

    allocate( start_positions(3, my_lattice%dim_lat(1), my_lattice%dim_lat(2), my_lattice%dim_lat(3), size(my_motif%atomic)) )
    call get_position( start_positions, my_lattice%dim_lat, my_lattice%areal, my_motif )
    pos=reshape( start_positions, (/3, N_cell/) )
    deallocate(start_positions)
    call calculate_distances(distances,pos,my_lattice%areal,my_lattice%dim_lat,my_lattice%boundary)
    deallocate(pos)

    call associate_pointer(all_mode,my_lattice)
    call set_E_bandstructure(my_lattice%dim_mode,distances)

    call calculate_dispersion(all_mode, dispersion, my_lattice%dim_mode, nb_kpoints, N_cell)

!allocate(input_energy(10*size(dispersion)), DOS(10*size(dispersion)))
call read_params_DOS('input')
call init_Evector_DOS(input_energy)
call compute_DOS(dispersion, input_energy, DOS, N_cell)
!do i=1, N_cell
!    write(*,*) 'all_mode(', i, ')%w(:) = ', all_mode(i)%w(:)
!    write(*,*) 'my_lattice%dim_mode = ', my_lattice%dim_mode
!    write(*,*) ''
!    write(*,*) ''
!enddo
call check_norm_wavefct(all_mode, my_lattice%dim_mode)

!io=open_file_write('dispersion_energy.txt')
!do i=1,nb_kpoints
!  write(io,'(2(f16.10,2x))') real(dispersion(i)),aimag(dispersion(i))
!enddo
!call close_file('dispersion_energy.txt',io)
!io=open_file_write('input_energy_DOS.txt')
!do i=1,size(input_energy)
!  write(io,'(2(f16.10,2x))') real(input_energy(i)),aimag(input_energy(i))
!enddo
!call close_file('input_energy_DOS.txt',io)
!io=open_file_write('DOS.txt')
!do i=1,size(DOS)
!  write(io,'(2(f16.10,2x))') real(DOS(i)),aimag(DOS(i))
!enddo
!call close_file('DOS.txt',io)
call rewrite_H_k(my_lattice%dim_mode)
!To get an easy processing with gnuplot, use
!paste kpoints data.txt | awk '{print $1 " " $4}' > dispersion.txt

end subroutine tightbinding