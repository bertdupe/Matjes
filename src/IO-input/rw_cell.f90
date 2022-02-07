module m_rw_cell
use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit, output_unit
implicit none
contains
subroutine read_cell(cell,lat)
    use m_derived_types
    use m_io_files_utils
    use m_io_utils
    use m_rw_motif
    type(lattice), intent(inout)    :: lat
    type(t_cell), intent(out)       :: cell
    ! internal
    type(atomtype),allocatable  ::  attype(:)
    integer :: io_input,n_mag,i
    
    io_input=open_file_read('input')
    call get_parameter(io_input,'input',attype)
    if(allocated(attype))then
        call get_parameter(io_input,'input',attype,cell)
        do i=1,size(cell%atomic)
            cell%atomic(i)%position=matmul(cell%atomic(i)%position,lat%areal)
        enddo
        call close_file('input',io_input)
        n_mag=count(cell%atomic(:)%moment.gt.0.0d0)

        write(output_unit,'(/A)') "Found the following atomic sites in the basic unit-cell:"
        write(output_unit,'(A)') "  name            position x (nm)   position y (nm)   position z (nm)       mag. moment      charge        mass     use phonon    number orbitals"
        do i=1,size(cell%atomic)
            write(output_unit,'(3X,A14,3(F12.8,6X),6X,3(F8.4,6X),L3,10X,I5)') cell%atomic(i)%name,cell%atomic(i)%position,cell%atomic(i)%moment,cell%atomic(i)%charge,cell%atomic(i)%mass,cell%atomic(i)%use_ph,cell%atomic(i)%orbitals
        enddo
        write(output_unit,'(/)')

    else 
        call close_file('input',io_input)
        WRITE(*,'(///A/A///)') "WARNING, FALLBACK TO OLD CELL INPUT VERSION","THIS MIGHT GIVE PROBLEMS WITH THE HAMILTONIAN DEFINITIONS"
        Call rw_motif(cell,lat)
    endif
end subroutine 


end module
