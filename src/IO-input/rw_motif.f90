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
    integer :: io_input,n_mag
    ! check the allocation of memory
    integer :: alloc_check,i,j,N
    integer :: natom,dimension(3)
    
    dimension=my_lattice%dim_lat

    
    io_input=open_file_read('input')
    call get_parameter(io_input,'input','motif',natom)
    allocate(my_motif%atomic(natom),stat=alloc_check)
    if (alloc_check.ne.0) write(6,'(a)') 'can not allocate motif'
    call get_parameter(io_input,'input','motif',natom,my_motif%atomic)
    call close_file('input',io_input)

    N=count(my_lattice%periodic)
    allocate(my_lattice%world(N),source=0)
    j=0
    do i=1,3
       if (my_lattice%periodic(N)) then
          j=j+1
          my_lattice%world(j)=my_lattice%dim_lat(i)
       endif
    enddo

    write(*,*) count(my_lattice%periodic)
    pause

end subroutine rw_motif


end module
