module m_rw_motif
implicit none
contains
subroutine rw_motif(my_motif,my_lattice)
    use m_derived_types
    use m_io_files_utils
    use m_io_utils
    type(lattice), intent(inout) :: my_lattice
    type(t_cell), intent(out) :: my_motif
    ! internal
    integer :: io_input,i,n_mag,i_at
    ! check the allocation of memory
    integer :: alloc_check
    integer :: natom,dimension(3)
    character(100)  ::  line
    integer         ::  stat
    
    dimension=my_lattice%dim_lat

    
    io_input=open_file_read('input')
    call get_parameter(io_input,'input','motif',natom)
    allocate(my_motif%atomic(natom),stat=alloc_check)
    if (alloc_check.ne.0) write(6,'(a)') 'can not allocate motif'
    call get_parameter(io_input,'input','motif',natom,my_motif%atomic)
    call close_file('input',io_input)
    
    n_mag=count(my_motif%atomic(:)%moment.gt.0.0d0)
    ! size of the world
    if ((dimension(3).eq.1).and.(dimension(2).eq.1)) then
        allocate(my_lattice%world(1))
        my_lattice%world(1)=dimension(1)
        my_lattice%n_system=1
        if (n_mag.gt.1) my_lattice%n_system=12
    elseif (dimension(3).eq.1) then
        allocate(my_lattice%world(2))
        my_lattice%world=(/dimension(1),dimension(2)/)
        my_lattice%n_system=2
        if (n_mag.gt.1) my_lattice%n_system=22
    elseif ((dimension(3).eq.1).and.(dimension(2).eq.1).and.(dimension(1).eq.1)) then
        write(6,*) "dimension of the problem not correct"
        stop
    else
        allocate(my_lattice%world(3))
        my_lattice%world=(/dimension(1),dimension(2),dimension(3)/)
        my_lattice%n_system=3
        if (n_mag.gt.1) my_lattice%n_system=32
    endif
end subroutine rw_motif


end module
