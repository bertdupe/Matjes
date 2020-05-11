module m_init_config
use m_derived_types
use m_init_spiral
use m_init_DW
use m_init_Sk
use m_init_Sklattice
use m_init_constant
use m_io_utils
use m_init_heavyside
use m_io_files_utils
use mtprng, only : mtprng_state
use m_init_random_config
use m_lattice, only : my_order_parameters
use m_convert

private
public :: init_config
contains

!!!!!!!!!!!
! find the different starting configurations for the order parameters
subroutine init_config(fname,my_lattice,my_motif,state,ext_param)
implicit none
character(len=*), intent(in) :: fname
type (lattice), intent(inout) :: my_lattice
type (cell), intent(in) :: my_motif
type (simulation_parameters), intent(in) :: ext_param
type (mtprng_state),intent(inout) :: state
! internal variables
integer :: io,N_mode,index_mode,i
character(len=30) :: seed_name,mode_name,configuration
integer :: nconfig

nconfig=0
seed_name='configuration_'
io=open_file_read(fname)
N_mode=size(my_order_parameters)

do i=1,N_mode

configuration='constant'
mode_name=convert(seed_name,my_order_parameters(i)%name)
call get_parameter(io,fname,mode_name,configuration)

  select case (adjustl(configuration))
    case('spiral')
      call init_spiral(io,fname,my_lattice,my_motif,my_order_parameters(i)%name,my_order_parameters(i)%start,my_order_parameters(i)%end)
    case('domainwall')
      call init_DW(my_lattice,my_motif,my_order_parameters(i)%start,my_order_parameters(i)%end)
    case('heavyside')
      call init_heavyside(my_lattice,my_motif,my_order_parameters(i)%start,my_order_parameters(i)%end)
    case('skyrmion')
      call init_spiral(io,fname,my_lattice,my_motif,my_order_parameters(i)%name,my_order_parameters(i)%start,my_order_parameters(i)%end)
      call init_Sk(io,fname,my_lattice,my_motif,my_order_parameters(i)%name,my_order_parameters(i)%start,my_order_parameters(i)%end)
    case('skyrmionla')
      call init_spiral(io,fname,my_lattice,my_motif,my_order_parameters(i)%name,my_order_parameters(i)%start,my_order_parameters(i)%end)
      call init_Sk_lattice(io,fname,my_lattice,my_motif,my_order_parameters(i)%name,my_order_parameters(i)%start,my_order_parameters(i)%end)
    case('random')
      write(6,'(a)') 'random configuration was choosen'
      call init_random_config(my_lattice,my_motif,state,my_order_parameters(i)%start,my_order_parameters(i)%end)
    case default
      call init_constant_config(my_lattice,my_order_parameters(i)%name,my_order_parameters(i)%start,my_order_parameters(i)%end,ext_param)
  end select

enddo
call close_file(fname,io)

end subroutine init_config





!!!!!!!!!!!!!!!!!!!!!!!!!
! function that finds the index of the mode of interest
!!!!!!!!!!!!!!!!!!!!!!!!!
integer function find_mode(mode_name)
implicit none
character(len=*), intent(in) :: mode_name
!internal
integer :: i,n_mode_total

find_mode=-1
n_mode_total=size(my_order_parameters)

do i=1,n_mode_total
   if (trim(my_order_parameters(i)%name).eq.trim(mode_name)) find_mode=i
enddo

if (find_mode.eq.-1) stop 'ERROR: Mode can not be found in initial-config'

end function find_mode

end module m_init_config
