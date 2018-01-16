module m_init_config
use m_derived_types
use m_init_spiral
use m_init_DW
use m_init_Sk
use m_init_Sklattice
use m_io_utils
use m_init_heavyside
implicit none
private
public :: init_config
contains

!!!!!!!!!!!
! find the different starting configurations for the order parameters
subroutine init_config(fname,my_lattice,my_motif)
character(len=*), intent(in) :: fname
type (lattice), intent(inout) :: my_lattice
type (cell), intent(in) :: my_motif
! internal variables
integer :: io
character(len=10) :: configuration
integer :: nconfig

nconfig=0
open(newunit=io,file=fname,form='formatted',status='old',action='read')

call get_parameter(io,fname,'configuration',configuration)

close(io)

select case (adjustl(configuration))
   case('spiral')
      call init_spiral(io,fname,my_lattice,my_motif)
   case('domainwall')
      call init_DW(my_lattice,my_motif)
   case('heavyside')
      call init_heavyside(my_lattice,my_motif)
   case('skyrmion')
      call init_Sk(io,fname,my_lattice,my_motif)
   case('skyrmionla')
      call init_Sk_lattice(io,fname,my_lattice,my_motif)
   case default
      write(6,'(a)') 'initial configuration not found'
      write(6,'(a)') 'the difference choices are'
      write(6,'(a)') 'spiral - domainwall - heavyside - skyrmion - skyrmionlattice'
      STOP
end select

end subroutine init_config

end module m_init_config
