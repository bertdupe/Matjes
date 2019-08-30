subroutine rw_motif(my_motif,my_lattice)
use m_derived_types
use m_io_files_utils
use m_io_utils
implicit none
! intent(inout)
type(lattice), intent(inout) :: my_lattice
type(cell), intent(out) :: my_motif
! internal
integer :: io_input,i
! check the allocation of memory
integer :: alloc_check
integer :: natom,dimension(3)

dimension=my_lattice%dim_lat

io_input=open_file_read('input')

call get_parameter(io_input,'input','motif',natom)

allocate(my_motif%atomic(natom),my_motif%i_mom(natom),stat=alloc_check)
if (alloc_check.ne.0) write(6,'(a)') 'can not allocate motif'

call get_parameter(io_input,'input','motif',natom,my_motif%atomic)

call close_file('input',io_input)

do i=1,natom
   if (abs(my_motif%atomic(i)%moment).lt.1.0d-8) then
      my_motif%i_mom(i)=.false.
   else
      my_motif%i_mom(i)=.true.
   endif
enddo

! size of the world
if ((dimension(3).eq.1).and.(dimension(2).eq.1)) then
    allocate(my_lattice%world(1))
    my_lattice%world(1)=dimension(1)
    my_lattice%n_system=1
    if (count(my_motif%i_mom).gt.1) my_lattice%n_system=12
elseif (dimension(3).eq.1) then
    allocate(my_lattice%world(2))
    my_lattice%world=(/dimension(1),dimension(2)/)
    my_lattice%n_system=2
    if (count(my_motif%i_mom).gt.1) my_lattice%n_system=22
elseif ((dimension(3).eq.1).and.(dimension(2).eq.1).and.(dimension(1).eq.1)) then
    write(6,*) "dimension of the problem not correct"
    stop
else
    allocate(my_lattice%world(3))
    my_lattice%world=(/dimension(1),dimension(2),dimension(3)/)
    my_lattice%n_system=3
    if (count(my_motif%i_mom).gt.1) my_lattice%n_system=32
endif


end subroutine rw_motif
