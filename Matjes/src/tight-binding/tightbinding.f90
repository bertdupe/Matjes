
subroutine tightbinding(lat,my_motif,io_simu,ext_param)
    use m_tb_params, only: TB_params, set_TB_params
    use m_basic_types, only : vec_point
    use m_derived_types, only : t_cell,lattice,io_parameter,simulation_parameters
    use m_operator_pointer_utils
    use m_lattice, only : my_order_parameters
    use m_tightbinding_r, only: tightbinding_r
!    use m_tightbinding_k, only: tightbinding_k

    implicit none
    ! internal parameter
    type(io_parameter), intent(in) :: io_simu
    type(lattice), intent(in) :: lat
    type(t_cell), intent(in) :: my_motif
    type(simulation_parameters), intent(in) :: ext_param

    ! N_cell is the variable that will contain the number of unit
    ! cells in the simulation
    integer :: N_cell

    integer :: i, TB_pos_ext(2)
    logical :: i_magnetic, i_TB
   
    N_cell=product(shape(lat%ordpar%l_modes))

    Call set_TB_params(N_cell)

    !do some initial testing real-space tight binding stuff
    if(TB_params%flow%do_r) Call tightbinding_r(lat,TB_params%H)   
   ! TAKE OUT K-SPACE STUFF SO FAR
   ! !do some initial testing reciprocal-space tight binding stuff
   ! if(TB_params%flow%do_k) Call tightbinding_k(TB_params%H,mode_magnetic,lat,my_motif)

end subroutine tightbinding
