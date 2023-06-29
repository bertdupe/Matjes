module m_energy_output_contribution
use m_H_public, only: t_H,energy_resolved
use m_hamiltonian_collection, only: hamiltonian
use m_derived_types, only: lattice
use m_convert, only: convert
use m_io_files_utils, only: open_file_write,close_file

private
public Eout_contrib_init,Eout_contrib_write
contains

subroutine Eout_contrib_init(H,io)
    type(hamiltonian), intent(in)   :: H
    integer, intent(out)            :: io

    integer     ::  i
    
    !write output file
    io=open_file_write('EM_energy_cont.dat')
    !!write header with information
    write(io,'(A)',advance="no") "# 1: step     2: time     3: total Energy     "
    do i=1,H%size_H()
        write(io,'(I4,A,A,5X)',advance="no") i+3,': ', trim(adjustl(H%get_desc(i)))
    enddo
    write(io,'(A)') " "
end subroutine

subroutine Eout_contrib_write(H,step,time,lat,io)
    type(hamiltonian),intent(inout) :: H
    integer,intent(in)              :: step
    real(8),intent(in)              :: time
    type(lattice),intent(inout)     :: lat
    integer, intent(in)             :: io

    real(8),allocatable :: E(:)
    integer             :: i, N_H

    N_H=H%size_H()
    allocate(E(N_H))
    Call H%energy_resolved(lat,E)
    write(*,*) 'in eout_contrib_write, E=',E(:), ' Etot=',sum(E)
    write(io,'(I6,2E16.8)',advance='no') step,time,sum(E) 
    do i=1,N_H
        write(io,'(E16.8)',advance='no') E(i) 
    end do
    write(io,'(A)') " " 
end subroutine
    
end module

