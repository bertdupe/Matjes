module m_excitation_from_file
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit
use m_io_utils
use m_io_files_utils

contains

subroutine read_excitation_from_file(fname,excitations_from_file)
  character(len=*),intent(in)       ::  fname
  real(8),allocatable, intent(out)  :: excitations_from_file(:,:)

  integer :: io_input,i,N_spin

  if (allocated(excitations_from_file)) then
     write(output_unit,'(a)') 'excitations_from_file already allocated'
     return
  endif

  N_spin=get_lines(fname)
  allocate(excitations_from_file(3,N_spin),source=0.0d0)

  io_input=open_file_read(fname)
  do i=1,N_spin
     read(io_input,*) excitations_from_file(:,i)
  enddo
  call close_file(fname,io_input)

end subroutine


end module m_excitation_from_file
