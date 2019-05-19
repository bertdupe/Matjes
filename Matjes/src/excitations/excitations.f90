module m_excitations
use m_basic_types, only : excitations
use m_simu_parameters, only : type_excitations
use m_rampe
use m_heavyside
use m_TPulse

type cycle
    type(excitations), allocatable :: external_param(:)
    integer :: counter
    character(len=30) :: name
end type cycle

type(cycle) :: EM_of_t

logical :: i_excitations

private
public :: get_excitations,update_EM_fields

contains

!
! routine that chooses if you are doing a EM_cycle, changing temperature or something else
!
subroutine get_excitations(fname)
use m_io_utils
use m_io_files_utils
implicit none
character(len=*), intent(in) :: fname
!internal variable
integer :: io_input,n_variable,i,n_excite
logical :: test
character(len=30) :: excitations

test=.false.
i_excitations=.false.
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

if (.not.i_excitations) then
   write(6,'(a)') 'no excitations found'
   call close_file('input',io_input)
   return
endif

EM_of_t%counter=0
EM_of_t%name=excitations

select case (excitations)
! part of the rampe function
  case('rampe')

    EM_of_t%counter=get_num_variable(io_input,EM_of_t%name)

    allocate(EM_of_t%external_param(EM_of_t%counter))

    call get_parameter(io_input,EM_of_t%name,trim(excitations),EM_of_t%external_param)

! part of the step function
  case('heavyside')

    EM_of_t%counter=get_num_variable(io_input,EM_of_t%name)

    allocate(EM_of_t%external_param(EM_of_t%counter))

    call get_parameter(io_input,EM_of_t%name,trim(excitations),EM_of_t%external_param)

  case('TPulse')

    write(6,'(a)') 'initialization of the electron temperature pulse'
    call get_parameter_TPulse(io_input,'input')

  case default
     write(6,'(a)') 'no excitations selected'

end select

call close_file('input',io_input)

end subroutine get_excitations

!
! routine that updates the external parameters
!
subroutine update_EM_fields(time,kt,h_int,check)
implicit none
integer, intent(in) :: time
real(kind=8), intent(inout) :: kt,h_int(:),check(:)
! internal parameters
integer :: size_excitations

if (.not.i_excitations) return

size_excitations=size(EM_of_t%external_param)

select case (EM_of_t%name)
  case('rampe')

    if (EM_of_t%counter.le.size_excitations) call update_rampe(time,kt,h_int,EM_of_t%external_param,EM_of_t%counter,check)

  case('heavyside')

    if (EM_of_t%counter.le.size_excitations) call update_heavyside(time,kt,h_int,EM_of_t%external_param,EM_of_t%counter,check)

  case('TPulse')

    if (EM_of_t%counter.le.size_excitations) call update_TPulse(time,kt)

  case default
    stop 'excitation not implemented'

end select

check=0.0d0

end subroutine update_EM_fields

!
! get the number of varibales on which you are doing the rampe
!

function get_num_variable(io,vname)
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
   if ( str(1:len_string) == vname) then
      backspace(io)
      read(io,*) dummy,dum_logic,test
      test=trim(adjustl(test))
      if ((test.eq.'H_ext').or.(test.eq.'E_ext').or.(test.eq.'temp')) then
         looping=.true.
         counter=counter+1
      endif
      do while (looping)
         looping=.false.
         read(io,*) test
         test=trim(adjustl(test))
         if ((test.eq.'H_ext').or.(test.eq.'E_ext').or.(test.eq.'temp')) then
            looping=.true.
            counter=counter+1
         endif
      enddo
   endif
enddo

get_num_variable=counter

end function get_num_variable

end module m_excitations
