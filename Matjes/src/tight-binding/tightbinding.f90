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
    use m_derived_types, only : cell,lattice,io_parameter,simulation_parameters

    implicit none
    ! internal parameter
    type(io_parameter), intent(in) :: io_simu
    type(lattice), intent(in) :: my_lattice
    type(cell), intent(in) :: my_motif
    type(simulation_parameters), intent(in) :: ext_param
end subroutine tightbinding
