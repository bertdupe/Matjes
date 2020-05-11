module m_excitations
use m_simu_parameters, only : type_excitations
use m_rampe
use m_heavyside
use m_TPulse
use m_EMwave
use m_shape_excitations
use m_get_position
use m_convert

! variable that contains the the excitations form (sweep of EM field...)
type excitations
    real(kind=8) :: start_value(3),end_value(3)
    real(kind=8) :: t_start,t_end
! name of the order parameter to change
    character(len=30) :: name
end type excitations

type excitation_order_parameter
! field concerned by the order parameter
  real(kind=8), allocatable :: field(:)
! name of the shape in field
  character(len=30) :: name
! pointer towards the norm of the field
  procedure(shape_norm), pointer, nopass :: norm => null()
end type excitation_order_parameter

type cycle_excitations
    type(excitations), allocatable :: temporal_param(:)
    type(excitation_order_parameter), allocatable :: spatial_param(:)
    ! count the excitations
    integer :: counter
    ! name of excitations
    character(len=30), allocatable :: name(:)
end type cycle_excitations





!
! position of the modes on the lattice
!
real(kind=8), allocatable :: position(:,:)
!
! EM cycle
!
type(cycle_excitations) :: EM_of_t

logical :: i_excitations

private
public :: get_excitations,update_ext_EM_fields,associate_excitation,update_EMT_of_r

contains

!
! routine that chooses if you are doing a EM_cycle, changing temperature or something else
!
subroutine get_excitations(fname,excitation,input_excitations)
use m_io_utils
use m_io_files_utils
use m_convert
implicit none
integer, intent(inout) :: input_excitations
character(len=*), intent(in) :: fname
logical, intent(out) :: excitation
!internal variable
integer :: io_input,i,n_excite,N,j
logical :: test
character(len=30) :: excitations
character(len=30) :: name_variable

test=.false.
i_excitations=.false.
excitation=.false.
input_excitations=0
n_excite=size(type_excitations)

io_input=open_file_read(fname)

call get_parameter(io_input,fname,'num_excitations',input_excitations)

allocate(EM_of_t%name(input_excitations),EM_of_t%temporal_param(input_excitations),EM_of_t%spatial_param(input_excitations))
EM_of_t%counter=input_excitations

if (input_excitations.ne.0) then
  i_excitations=.true.
  excitation=.true.
endif

if (.not.i_excitations) then
   write(6,'(a)') 'no excitations found'
   call close_file(fname,io_input)
   return
endif

do i=1,EM_of_t%counter
  do j=1,n_excite
    test=.false.
    name_variable=convert(trim(type_excitations(j)%name),'_')
    name_variable=convert(trim(name_variable),i)
    call get_parameter(io_input,fname,name_variable,test)
    if (test) EM_of_t%name(i)=type_excitations(j)%name
  enddo
enddo

do j=1,EM_of_t%counter
 excitations=EM_of_t%name(j)

 select case (excitations)
! part of the rampe function
  case('rampe')
    call read_excitations(io_input,fname,excitations,EM_of_t%temporal_param(j))

! part of the step function
  case('heavyside')

     call read_excitations(io_input,fname,excitations,EM_of_t%temporal_param(j))

  case('TPulse')

    write(6,'(a)') 'initialization of the electron temperature pulse'
    call get_parameter_TPulse(io_input,fname)

  case('EMwave')

    write(6,'(a)') 'initialization of the electron temperature pulse'
    call get_parameter_EMwave(io_input,fname)

  case default
     write(6,'(a)') 'no excitations selected'

 end select

 call get_shape(io_input,fname,EM_of_t%spatial_param(j)%name,EM_of_t%spatial_param(j)%norm)

 if ((EM_of_t%temporal_param(j)%name.eq.'Bfield').or.(EM_of_t%temporal_param(j)%name.eq.'Efield')) then
   allocate(EM_of_t%spatial_param(j)%field(3))
 else
   allocate(EM_of_t%spatial_param(j)%field(1))
 endif
 EM_of_t%spatial_param(j)%field=0.0d0
enddo

call close_file(fname,io_input)

N=get_lines('positions.dat')
allocate(position(3,N))
position=0.0d0
call get_position(position,'positions.dat')

!
end subroutine get_excitations

!
! routine that associates the field to the table of pointer that should be used
!
subroutine associate_excitation(i_excite,point,all_mode,my_order_parameters)
use m_derived_types, only : order_parameter
use m_basic_types, only : vec_point
implicit none
type(order_parameter), intent(in) :: my_order_parameters(:)
type(vec_point), intent(in):: all_mode(:)
type(vec_point), intent(inout):: point(:)
integer, intent(in) :: i_excite
! internal variables
integer :: i,j

do i=1,size(my_order_parameters)
   if( my_order_parameters(i)%name.eq.EM_of_t%temporal_param(i_excite)%name) then
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
type(vec_point), intent(inout) :: mode(:,:)
! internal
real(kind=8) :: r(3)
integer :: j

r=position(:,i)

do j=1,EM_of_t%counter
  mode(i,j)%w=EM_of_t%spatial_param(j)%field*EM_of_t%spatial_param(j)%norm(r,shape_excitation%center,shape_excitation%cutoff)
enddo

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
integer :: size_excitations,size_field,j
character(len=30) :: form
real(kind=8),allocatable :: field_ini(:)

if (.not.i_excitations) return

size_excitations=size(EM_of_t%temporal_param)

do j=1,EM_of_t%counter

  size_field=size(EM_of_t%spatial_param(j)%field)
  allocate(field_ini(size_field))
  field_ini=EM_of_t%spatial_param(j)%field

  select case (EM_of_t%name(j))
    case('rampe')

      if (EM_of_t%counter.le.size_excitations) call update_rampe(time,EM_of_t%spatial_param(j)%field, &
      & EM_of_t%temporal_param(j)%t_start,EM_of_t%temporal_param(j)%t_end,EM_of_t%temporal_param(j)%start_value, &
      & EM_of_t%temporal_param(j)%end_value,EM_of_t%counter)

    case('heavyside')

      if (EM_of_t%counter.le.size_excitations) call update_heavyside(time,EM_of_t%spatial_param(j)%field, &
      & EM_of_t%temporal_param(j)%t_start,EM_of_t%temporal_param(j)%t_end,EM_of_t%temporal_param(j)%start_value, &
      & EM_of_t%temporal_param(j)%end_value,EM_of_t%counter)

    case('TPulse')

      if (EM_of_t%counter.le.size_excitations) call update_TPulse(time,EM_of_t%spatial_param(j)%field)

    case('EMwave')

      if (EM_of_t%counter.le.size_excitations) call update_EMwave(time,EM_of_t%spatial_param(j)%field)

    case default
      stop 'excitation not implemented'

  end select

  if (norm(field_ini-EM_of_t%spatial_param(j)%field).gt.1.0d-5) then
    form=convert('(a,',size_field,'f14.6/)')
    write(6,'(/a,2x,f16.6)') 'real time',time
    write(6,form) 'field value ',EM_of_t%spatial_param(j)%field(:)
  endif
  if (abs(check(2)).gt.1.0d-8) write(6,'(a,2x,f16.6)') 'Final Temp', check(1)/check(2)/2.0d0/k_B

  deallocate(field_ini)

enddo

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

!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get the excitations for all the cycles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
subroutine read_excitations(io,fname,vname,excite)
use m_io_utils, only : check_read
implicit none
type(excitations), intent(inout) :: excite
integer, intent(in) :: io
character(len=*), intent(in) :: vname,fname
! internal variable
integer :: fin,len_string,nread,check
character(len=100) :: str
character(len=100) :: dummy
logical :: dum_logic

nread=0
len_string=len(trim(adjustl(vname)))


excite%start_value(:)=0.0d0
excite%end_value(:)=0.0d0
excite%t_start=0
excite%t_end=0
excite%name=fname


rewind(io)
do
   read (io,'(a)',iostat=fin) str
   if (fin /= 0) exit
   str= trim(adjustl(str))

   if (len_trim(str)==0) cycle
   if (str(1:1) == '#' ) cycle

!cccc We start to read the input
   if ( str(1:len_string) == trim(adjustl(vname)) ) then
      nread=nread+1
      backspace(io)
      read(io,*) dummy,dum_logic,excite%name
      if ((excite%name.eq.'Bfield').or.(excite%name.eq.'Efield')) then
         backspace(io)
         read(io,*) dummy,dum_logic,dummy,excite%start_value(1:3),excite%end_value(1:3),excite%t_start,excite%t_end
      else
         backspace(io)
         read(io,*) dummy,dum_logic,dummy,excite%start_value(1),excite%end_value(1),excite%t_start,excite%t_end
      endif

   endif

enddo

check=check_read(nread,vname,fname)

end subroutine read_excitations

end module m_excitations
