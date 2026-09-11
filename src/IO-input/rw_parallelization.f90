module m_parallelization_io
use mpi_basic,only: mpi_type,mpi_distv
use m_io_utils
use m_io_files_utils
use mpi_util
use, intrinsic  ::  ISO_FORTRAN_ENV, only: output_unit
implicit none


contains






subroutine rw_param_paral(mpi_outer,mpi_inner,user_def,mpi_in)
  type(mpi_distv),intent(out)     :: mpi_outer
  type(mpi_type) ,intent(out)     :: mpi_inner
  logical        ,intent(inout)   :: user_def
  type(mpi_type) ,intent(in)      :: mpi_in

  integer :: io_input

  if (mpi_in%ismas) then

     io_input=open_file_read('input')

     call get_parameter(io_input,'input','Np_outer',mpi_outer%Np)
! the Hamiltonian is parallelized on the inner processes
     call get_parameter(io_input,'input','Np_inner',mpi_inner%Np)

     call close_file('input',io_input)

     if (mpi_outer%Np*mpi_inner%Np.ne.1) then
        write(6,'(/a)') 'user defined parallelization'
        write(6,'(a,I4,a/)') 'The Hamiltonian will be parallelized on ', mpi_inner%Np,'  procs'
        user_def=.true.
     else
        write(6,'(/a/)') 'automatic parallelization'
     endif
  endif

  call bcast(mpi_outer%Np,mpi_in)
  call bcast(mpi_inner%Np,mpi_in)
  call bcast(user_def,mpi_in)

end subroutine

end module
