module m_energy_output_contribution
use m_H_public, only: t_H,energy_resolved
use m_derived_types, only: lattice
use m_convert, only: convert
use m_io_files_utils, only: open_file_write,close_file

private
public Eout_contrib_init,Eout_contrib_write
contains

subroutine Eout_contrib_init(Hams,io)
    class(t_H), intent(in)          :: Hams(:)
    integer, intent(out)            :: io

    integer     ::  i
    
    !write output file
    io=open_file_write('EM_energy_cont.dat')
    !!write header with information
    write(io,'(A)',advance="no") "# 1: step     2: time     3: Energy average     "
    do i=1,size(Hams)
        write(io,'(I4,A,A,5X)',advance="no") i+3,': ', trim(adjustl(Hams(i)%desc))
    enddo
    write(io,'(A)') " "
end subroutine

subroutine Eout_contrib_write(step,time,Hams,lat,io)
    integer,intent(in)              :: step
    real(8),intent(in)              :: time
    class(t_H), intent(in)          :: Hams(:)
    class(lattice),intent(in)       :: lat
    integer, intent(in)             :: io

    real(8)         :: E(size(Hams))
    integer         :: i

    Call energy_resolved(hams,lat,E)
    write(io,'(I6,2E16.8)',advance='no') step,time,sum(E) 
    do i=1,size(E)
        write(io,'(E16.8)',advance='no') E(i) 
    end do
    write(io,'(A)') " " 
end subroutine
    
end module

