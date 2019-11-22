module m_excitations
use m_basic_types, only : excitations
use m_simu_parameters, only : type_excitations
use m_rampe
use m_heavyside
use m_TPulse
use m_EMwave
use m_shape_excitations
use m_get_position
use m_convert

type cycle
    type(excitations), allocatable :: external_param(:)
    real(kind=8), allocatable :: field(:)
    integer :: counter
    character(len=30) :: name
    procedure(shape_norm), pointer, nopass :: norm => null()
end type cycle

real(kind=8), allocatable :: position(:,:)

type(cycle) :: EM_of_t

logical :: i_excitations

private
public :: get_excitations,update_ext_EM_fields,associate_excitation,update_EMT_of_r

contains

!
! routine that chooses if you are doing a EM_cycle, changing temperature or something else
!
subroutine get_excitations(fname,excitation)
use m_io_utils
use m_io_files_utils
use m_convert
implicit none
character(len=*), intent(in) :: fname
logical, intent(out) :: excitation
!internal variable
integer :: io_input,n_variable,i,n_excite,N
logical :: test
character(len=30) :: excitations
character(len=30) :: name_variable

test=.false.
i_excitations=.false.
excitation=.false.
n_excite=size(type_excitations)

io_input=open_file_read(fname)

do i=1,n_excite
   call get_parameter(io_input,fname,trim(type_excitations(i)%name),test)
   if (test) then
      excitations=trim(type_excitations(i)%name)
      exit
   endif
enddo

i_excitations=test
excitation=test

if (.not.i_excitations) then
   write(6,'(a)') 'no excitations found'
   call close_file(fname,io_input)
   return
endif

EM_of_t%counter=0
EM_of_t%name=excitations

select case (excitations)
! part of the rampe function
  case('rampe')

    EM_of_t%counter=get_num_variable(io_input,EM_of_t%name)

    allocate(EM_of_t%external_param(EM_of_t%counter))

    do i=1, EM_of_t%counter
      name_variable=convert(trim(excitations),'_',EM_of_t%name,'_',i)
      call get_parameter(io_input,EM_of_t%name,name_variable,EM_of_t%external_param(i))
    enddo

! part of the step function
  case('heavyside')

    EM_of_t%counter=get_num_variable(io_input,EM_of_t%name)

    allocate(EM_of_t%external_param(EM_of_t%counter))

    do i=1, EM_of_t%counter
      name_variable=convert(EM_of_t%name,'_')
      excitations=convert(name_variable,i)
      call get_parameter(io_input,fname,excitations,EM_of_t%external_param(i))
    enddo

  case('TPulse')

    write(6,'(a)') 'initialization of the electron temperature pulse'
    call get_parameter_TPulse(io_input,fname)

  case('EMwave')

    write(6,'(a)') 'initialization of the electron temperature pulse'
    call get_parameter_EMwave(io_input,fname)

  case default
     write(6,'(a)') 'no excitations selected'

end select

call get_shape(io_input,fname,EM_of_t%external_param(1)%name,EM_of_t%norm)

call close_file(fname,io_input)

if ((EM_of_t%external_param(1)%name.eq.'Hfield').or.(EM_of_t%external_param(1)%name.eq.'Efield')) then
   allocate(EM_of_t%field(3))
else
   allocate(EM_of_t%field(1))
endif

EM_of_t%field=0.0d0

N=get_lines('positions.dat')
allocate(position(3,N))
position=0.0d0
call get_position(position,'positions.dat')

end subroutine get_excitations

!
! routine that associates the field to the table of pointer that should be used
!
subroutine associate_excitation(point,all_mode,my_order_parameters)
use m_derived_types, only : order_parameter
use m_basic_types, only : vec_point
implicit none
type(order_parameter), intent(in) :: my_order_parameters(:)
type(vec_point), intent(in):: all_mode(:)
type(vec_point), intent(inout):: point(:)
! internal variables
integer :: i,j

do i=1,size(my_order_parameters)
   if( my_order_parameters(i)%name.eq.EM_of_t%external_param(1)%name) then
      write(6,'(2a)') 'excitations was found on  ', my_order_parameters(i)%name
      do j=1,size(all_mode)
         point(j)%w => all_mode(j)%w(my_order_parameters(i)%start:my_order_parameters(i)%end)
      enddo
   endif
enddo

write(6,'(/,a,/)') 'the pointers for the excitations has been allocated in associate_excitation'

end subroutine

!
! routine that updates the the field depending on the position
!
subroutine update_EMT_of_r(i,mode)
use m_basic_types, only : vec_point
implicit none
integer, intent(in) :: i
real(kind=8), intent(inout) :: mode(3)
! internal
real(kind=8) :: r(3)

r=position(:,i)

mode=EM_of_t%field*EM_of_t%norm(r,shape_excitation%center,shape_excitation%cutoff)

end subroutine





!
! routine that updates the external parameters
!
subroutine update_ext_EM_fields(time,check)
use m_constants, only : k_b
use m_vector, only : norm
implicit none
real(kind=8), intent(in) :: time
real(kind=8), intent(inout) :: check(:)
! internal parameters
integer :: size_excitations,size_field
character(len=30) :: form
real(kind=8),allocatable :: field_ini(:)

if (.not.i_excitations) return

size_excitations=size(EM_of_t%external_param)
size_field=size(EM_of_t%field)
allocate(field_ini(size_field))
field_ini=EM_of_t%field

select case (EM_of_t%name)
  case('rampe')

    if (EM_of_t%counter.le.size_excitations) call update_rampe(time,EM_of_t%field,EM_of_t%external_param,EM_of_t%counter)

  case('heavyside')

    if (EM_of_t%counter.le.size_excitations) call update_heavyside(time,EM_of_t%field,EM_of_t%external_param,EM_of_t%counter)

  case('TPulse')

    if (EM_of_t%counter.le.size_excitations) call update_TPulse(time,EM_of_t%field)

  case('EMwave')

    if (EM_of_t%counter.le.size_excitations) call update_EMwave(time,EM_of_t%field)

  case default
    stop 'excitation not implemented'

end select

if (norm(field_ini-EM_of_t%field).gt.1.0d-5) then
   form=convert('(a,',size_field,'f14.6/)')
   write(6,'(/a,2x,f16.6)') 'real time',time
   write(6,form) 'field value ',EM_of_t%field(:)
endif
if (abs(check(2)).gt.1.0d-8) write(6,'(a,2x,f16.6)') 'Final Temp', check(1)/check(2)/2.0d0/k_B

check=0.0d0

end subroutine update_ext_EM_fields

!
! get the number of varibales on which you are doing the rampe
!

function get_num_variable(io,vname)
use m_io_utils, only : check_last_char
implicit none
integer, intent(in) :: io
character(len=*) :: vname
integer :: get_num_variable
! internal variable
integer :: counter,len_string,fin
character(len=100) :: str,test,dummy
logical :: looping,dum_logic

looping=.false.
counter=0
len_string=len(trim(adjustl(vname)))

rewind(io)
do
   read (io,'(a)',iostat=fin) str
   if (fin /= 0) exit
   str= trim(adjustl(str))

   if (len_trim(str)==0) cycle
   if (str(1:1) == '#' ) cycle

!cccc We start to read the input
   if (str(1:len_string) == vname)  then
      backspace(io)
      read(io,*) dummy,dum_logic,test
      test=trim(adjustl(test))
      if ((test.eq.'Bfield').or.(test.eq.'Efield').or.(test.eq.'temp')) then
         looping=.true.
         counter=counter+1
      endif
   endif
enddo

if (counter.eq.0) stop 'get_num_variables could not read Bfield or Efield or temp'

get_num_variable=counter

end function get_num_variable

end module m_excitations
