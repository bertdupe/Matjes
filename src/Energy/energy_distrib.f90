module m_energyfield
use m_hamiltonian_collection, only: hamiltonian
use m_derived_types, only: lattice
use m_convert, only: convert
use m_io_files_utils, only: open_file_write,close_file
use,intrinsic :: iso_c_binding, only: c_bool
private
public :: write_energy_field

contains

subroutine write_energy_field(tag,H,lat,order,ismas)
    !write the cell-resolved energy terms of Hams out to external file
    integer, intent(in)             :: tag
    type(hamiltonian),intent(inout) :: H
    type(lattice), intent(inout)    :: lat
    integer,intent(in)              :: order
    logical(c_bool),intent(in)      :: ismas
    !internal
    real(8),allocatable     :: energies(:,:)
    character(len=30)       :: fname,rw_format
    integer                 :: i,io, NH

    Call H%energy_distrib(lat,order,energies)
    if(ismas)then
        !write output file
        NH=H%size_H()
        write(rw_format,'(a,I4,a)') '(',NH,'E16.8)'
        fname=convert('EnDistrib_',tag,'.dat')
        io=open_file_write(fname)
        !!write header with information
        write(io,'(A)',advance="no") "# "
        do i=1,NH
            write(io,'(I4,A,A,5X)',advance="no") i,': ', trim(adjustl(H%get_desc(i)))
        enddo
        write(io,'(A)') " "
        !!write data
        write(io,rw_format) transpose(energies)
        call close_file(fname,io)
    endif
end subroutine

end module
