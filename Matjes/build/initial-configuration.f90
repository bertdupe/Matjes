module m_init_config
use m_derived_types
use m_init_spiral
use m_init_DW
use m_init_Sk
use m_init_Sklattice
use m_io_utils
use m_init_heavyside
use m_io_files_utils
use mtprng, only : mtprng_state
use m_init_random_config

private
public :: init_config
contains

!!!!!!!!!!!
! find the different starting configurations for the order parameters
subroutine init_config(fname,my_lattice,my_motif,state)
implicit none
character(len=*), intent(in) :: fname
type (lattice), intent(inout) :: my_lattice
type (cell), intent(in) :: my_motif
type (mtprng_state),intent(inout) :: state
! internal variables
integer :: io
character(len=10) :: configuration
integer :: nconfig

nconfig=0
io=open_file_read(fname)

call get_parameter(io,fname,'configuration',configuration)

select case (adjustl(configuration))
   case('spiral')
      call init_spiral(io,fname,my_lattice,my_motif)
   case('domainwall')
      call init_DW(my_lattice,my_motif)
   case('heavyside')
      call init_heavyside(my_lattice,my_motif)
   case('skyrmion')
      call init_spiral(io,fname,my_lattice,my_motif)
      call init_Sk(io,fname,my_lattice,my_motif)
   case('skyrmionla')
      call init_spiral(io,fname,my_lattice,my_motif)
      call init_Sk_lattice(io,fname,my_lattice,my_motif)
   case default
      write(6,'(a)') 'random configuration was choosen'
      call init_random_config(my_lattice,my_motif,state)
end select

call close_file(fname,io)

end subroutine init_config

end module m_init_config
