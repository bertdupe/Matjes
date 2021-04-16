module m_energyfield
use m_H_public, only: t_H
use m_derived_types, only: lattice
use m_convert, only: convert
use m_io_files_utils, only: open_file_write,close_file
private
public :: write_energy_field

contains

subroutine write_energy_field(tag,Hams,lat,order)
    !write the cell-resolved energy terms of Hams out to external file
    integer, intent(in)             :: tag
    class(t_H), intent(in)          :: Hams(:)
    type(lattice), intent(inout)    :: lat
    integer,intent(in)              :: order
    !internal
    real(8),allocatable     :: energies(:,:)
    character(len=30)       :: fname,rw_format
    integer                 :: i,io
    !calculate the energy values
    allocate(energies(lat%Ncell*lat%site_per_cell(order),size(Hams)),source=0.0d0)
    do i=1,size(Hams)
        Call hams(i)%energy_dist(lat,order,energies(:,i))
    enddo
    !write output file
    write(rw_format,'(a,I4,a)') '(',size(hams),'E16.8)'
    fname=convert('EnDistrib_',tag,'.dat')
    io=open_file_write(fname)
    !!write header with information
    write(io,'(A)',advance="no") "# "
    do i=1,size(Hams)
        write(io,'(I4,A,A,5X)',advance="no") i,': ', trim(adjustl(Hams(i)%desc))
    enddo
    write(io,'(A)') " "
    !!write data
    write(io,rw_format) transpose(energies)
    call close_file(fname,io)
end subroutine

end module
