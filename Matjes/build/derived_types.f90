       module m_derived_types

! type vector
       type vec
          real(kind=8) :: w(3)
       end type vec

       type vec_point
          real(kind=8), pointer :: w(:)
       end type vec_point


! Order parameter of the simulations
       type Xvec
        type(vec), allocatable, dimension(:) :: X
       end type Xvec

! this target is for the boolean variable type
type bool_var
     logical :: value
     character(len=:), allocatable :: name
     character(len=:), allocatable :: var_name
end type bool_var


! unit cells
! the unit cell can be magnetic, ferroelectric or nothing and can also have transport
type cell
     real(kind=8), allocatable :: mom(:), pos(:,:)
     logical, allocatable :: i_mom(:)
end type cell

type lattice
     real(kind=8) :: areal(3,3),astar(3,3),alat(3)
     integer :: dim_lat(3),n_system
     integer, allocatable :: world(:)
     logical :: boundary(3)
     real(kind=8), pointer :: l_modes(:,:,:,:)
end type lattice

! parameters for printing in and out
       type io_parameter
        logical :: io_dispersion,io_qorien
! Go in the spmstm program of Tobias
        Logical :: io_spstmL, io_spstmonly
! plot the spin structure (or the order parameter structure)
        Logical :: io_Xstruct
! plot the fourrier tranform of the spin structure (or the order parameter structure)
        Logical :: io_fft_Xstruct
! plot the topological charge density distribution
        Logical :: io_topo
! plot the emergent magnetic field
        logical :: i_topohall
! frequency of the plotting
        integer :: io_frequency
! plot the warnings or not
        logical :: io_warning
! frequency of writting of the data
        integer :: io_writing
       end type io_parameter






! parameters for the type of simulations that you are running
       type type_simu
        logical :: i_pimc,i_dynamic,i_metropolis,i_gneb,i_paratemp,i_minimization,i_entropic,i_r_texture,i_TB
       end type type_simu




!Operators of the simulations
       type Operator
        class (*), pointer, dimension(:,:) :: Op
        integer :: il,ic
        integer :: nline,ncolumn
       end type Operator


!!!!
!Hamiltonian of the simulations
!!!!

       type H_loc
        real(kind=8), pointer :: H_loc(:,:)
       end type H_loc

       type Hamiltonian
        class (*), pointer, dimension(:) :: Ham
        integer, dimension(:), allocatable :: Op_order
        integer :: N_op
       end type Hamiltonian

!Hamiltonian of the simulations
       type real_local_Hamil
        real(kind=16), allocatable, dimension(:,:) :: coeff
       end type real_local_Hamil

       type im_local_Hamil
        complex(kind=16), allocatable, dimension(:,:) :: coeff
       end type im_local_Hamil


       end module m_derived_types
