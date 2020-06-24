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
    use m_lattice, only : my_order_parameters

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
    type(vec_point),allocatable :: all_mode(:), mode_magnetic(:), mode_TB(:)
    ! sense gives the sense of the FFT transform
    real(kind=8) :: sense
    ! array containing all the positions in the lattice
    real(kind=8), allocatable :: start_positions(:,:,:,:,:),pos(:,:),distances(:,:)
    ! Etot gives the total energy contained in the system, kt is the thermal energy
    real(kind=8) :: Etot, kt
    ! E_F gives the Fermi energy
    real(kind=8) :: E_F
    ! eps_nk is a vector containing all the eigenvalues
    complex(kind=16), allocatable :: eps_nk(:)
    real(kind=8), allocatable :: eigval(:,:)

    complex(kind=16), allocatable :: dispersion(:), input_energy(:)
    integer :: io, i, nb_kpoints, TB_pos_start, TB_pos_end
    real(kind=8) :: N_electrons
    logical :: i_magnetic, i_TB

    kt=0.0d0
!    N_electrons=0.0d0

    call read_params_DOS('input')

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

! Allocating the different variables
!
    call associate_pointer(all_mode,my_lattice)
! magnetization
do i=1,size(my_order_parameters)
  if ('magnetic'.eq.trim(my_order_parameters(i)%name)) then
   allocate(mode_magnetic(N_cell))
   call dissociate(mode_magnetic,N_cell)
   call associate_pointer(mode_magnetic,all_mode,'magnetic',i_magnetic)
  endif

  if ('Tight-binding'.eq.trim(my_order_parameters(i)%name)) then
   allocate(mode_TB(N_cell))
   call dissociate(mode_TB,N_cell)
   call associate_pointer(mode_TB,all_mode,'Tight-binding',i_TB)

   TB_pos_start=my_order_parameters(i)%start
   TB_pos_end=my_order_parameters(i)%end
    allocate( eigval(N_cell*(TB_pos_end-TB_pos_start+1),nb_kpoints) )
    eigval = 0.0d0
  endif
enddo

!   initializing band structure and H(k)
call set_E_bandstructure(my_lattice%dim_mode,distances)
call rewrite_H_k(size(mode_TB(1)%w),TB_pos_start,TB_pos_end,distances)
deallocate(distances)


!The function "diagonalise_H_k" in file energy_k.f90 diagonalises the Hamiltonian
!for a given k-vector (it calls the function "Fourier_transform_H" inside)
!==> we need to loop over all the k-vectors to have the eigenenergies for
!all wavevectors
! les états vides ET les états pleins sont pris en compte
     do i=1,nb_kpoints
        call diagonalise_H_k(i, size(mode_TB(1)%w), -1.0d0, eigval(:,i))
     enddo


    N_electrons = check_norm_wavefct(all_mode, TB_pos_start, TB_pos_end, N_electrons)
    write(6,'(a,2x,f10.4)') ' N_electrons = ', N_electrons

!    do i=1,N_cell
!       E_F = 0.0d0
!       call compute_Fermi_level(eigval(i,:), N_electrons, E_F, kt)
!    enddo
!    call print_band_struct('N_bands.dat',eigval)
stop 'toto'
! diagonlisation uniquement avec les états
    call calculate_dispersion(all_mode, dispersion, my_lattice%dim_mode, nb_kpoints, N_cell)
    call print_band_struct('bands.dat',dispersion)

! faire la DOS
!allocate(input_energy(10*size(dispersion)), DOS(10*size(dispersion)))
    call init_Evector_DOS()
    call compute_DOS(dispersion, N_cell)
    call print_DOS('DOS.dat')



end subroutine tightbinding
