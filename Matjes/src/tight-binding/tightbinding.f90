
subroutine tightbinding(my_lattice,my_motif,io_simu,ext_param)
    use m_basic_types, only : vec_point
    use m_derived_types, only : cell,lattice,io_parameter,simulation_parameters
    use m_operator_pointer_utils
    use m_lattice, only : my_order_parameters
    use m_tightbinding_r, only: tightbinding_r
    use m_tightbinding_k, only: tightbinding_k
    use m_io_files_utils, only: open_file_read,close_file
    use m_io_utils, only: get_parameter

    implicit none
    ! internal parameter
    type(io_parameter), intent(in) :: io_simu
    type(lattice), intent(in) :: my_lattice
    type(cell), intent(in) :: my_motif
    type(simulation_parameters), intent(in) :: ext_param

    ! N_cell is the variable that will contain the number of unit
    ! cells in the simulation
    integer :: N_cell
    ! all_mode is an array of custom type vec_point that will
    ! contain all the modes present in the simulation
    type(vec_point),allocatable :: all_mode(:), mode_magnetic(:), mode_TB(:)

    integer :: i, TB_pos_ext(2),io_input
    logical :: i_magnetic, i_TB
    integer :: dimH
    logical :: do_TB_r,do_TB_k
    
    !get magnetization everywhere and set some pointers and control integers
    N_cell=product(shape(my_lattice%l_modes))
    allocate( all_mode(N_cell))
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
       TB_pos_ext(1)=my_order_parameters(i)%start
       TB_pos_ext(2)=my_order_parameters(i)%end
       dimH=N_cell*(TB_pos_ext(2)-TB_pos_ext(1)+1)
      endif
    enddo
#if CPP_SC
    if(.true.) dimH=dimH*2  !SC
#endif

    io_input=open_file_read('input')
    do_TB_r=.False.
    do_TB_k=.False.
    call get_parameter(io_input,'input','do_TB_k',do_TB_k)
    call get_parameter(io_input,'input','do_TB_r',do_TB_r)
    call close_file('input',io_input)
    !do some initial testing real-space tight binding stuff
    if(do_TB_r) Call tightbinding_r(dimH,TB_pos_ext,mode_magnetic)   
    !do some initial testing reciprocal-space tight binding stuff
    if(do_TB_k) Call tightbinding_k(dimH,TB_pos_ext,mode_magnetic,my_lattice,my_motif)

end subroutine tightbinding
