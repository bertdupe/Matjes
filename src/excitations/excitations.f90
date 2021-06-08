module m_excitations
use m_simu_parameters, only : type_excitations
use m_rampe
use m_heaviside
use m_TPulse
use m_EMwave
use m_shape_excitations
use m_get_position
use m_convert
use m_type_lattice
use m_excitation_norm
use m_excitation_shape
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit
implicit none

! variable that contains the the excitations form (sweep of EM field...)
type excitations
    real(8)           :: start_value(maxval(dim_modes_inner))
    real(8)           :: end_value(maxval(dim_modes_inner))
    integer           :: size_value
    real(kind=8)      :: t_start,t_end
! name of the order parameter to change
    character(len=30) :: name
    integer           :: op !order parameter as ordered by m_type_lattice
end type excitations

type excitation_order_parameter
! field concerned by the order parameter
  real(kind=8), allocatable :: field(:)
! name of the shape in field
  character(len=30)         :: name
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


type excitation
    integer     ::  op=0
    integer     ::  int_shape=0
    integer     ::  int_norm=0
end type

type excitation_combined
    type(excitation),allocatable        :: exc(:)
    type(excitation_norm),allocatable   :: norms(:)
    type(excitation_shape),allocatable  :: shapes(:)

    real(8), allocatable :: position(:,:)
end type

real(kind=8), allocatable :: position(:,:)
!
! position of the modes on the lattice
!
!
! EM cycle
!
type(cycle_excitations) :: EM_of_t

logical :: i_excitations

private
public :: set_excitations,update_ext_EM_fields,update_EMT_of_r
public :: read_all_excitations, update_exc, excitation_combined

contains

subroutine read_all_excitations(fname,excitations)
    use m_io_utils
    use m_io_files_utils
    use m_io_read_util
    use m_excitation_io
    character(len=*),intent(in)             :: fname
    type(excitation_combined),intent(inout) :: excitations

    type(excitation_io),allocatable     :: io_exc(:)

    integer             :: i,j, Nexc, io
   
    open(newunit=io,file=fname)
    !get different excitation norms
    Call read_excitation_norms(io,fname,excitations%norms)
    
    !get different excitation shapes
    Call read_excitation_shape(io,fname,excitations%shapes)

    !read excitation combinations
    Call read_excitation_io(io,fname,io_exc)

    Nexc=0
    if(allocated(io_exc)) Nexc=size(io_exc) 
    if(Nexc==0)then
        !no excitations inserted
        !clean up and leave
        if(allocated(excitations%shapes)) deallocate(excitations%shapes)
        if(allocated(excitations%norms )) deallocate(excitations%norms )
        return  
    else
        !check allocations
        if(.not.allocated(excitations%shapes))then
            write(error_unit,'(A)') "ERROR, Found some excitations, but 'excitation_shape' is not set"
            write(error_unit,'(A)') "CHECK INPUT"
            STOP
        endif
        if(.not.allocated(excitations%norms))then
            write(error_unit,'(A)') "ERROR, Found some excitations, but 'excitation_norm' is not set"
            write(error_unit,'(A)') "CHECK INPUT"
            STOP
        endif
    endif

    !set excitation entries
    allocate(excitations%exc(Nexc))
    do i=1,Nexc
        excitations%exc(i)%op=io_exc(i)%op
        !find associated shape entry
        do j=1,size(excitations%shapes)
            if(io_exc(i)%shape_name==excitations%shapes(j)%name)then
                if(excitations%shapes(j)%dim_mode/=dim_modes_inner(io_exc(i)%op))then
                    write(error_unit,'(A)') "ERROR, Excitation operator dimension does not fit associated excitation shape"
                    write(error_unit,'(2(AI3))') "Error associating excitation ",i,"with shape number",j
                    STOP
                endif
                excitations%exc(i)%int_shape=j
                exit
            endif
        enddo
        if(excitations%exc(i)%int_shape==0)then
            write(error_unit,'(A)') "ERROR, Excitation operator shape name not found in exitation shapes"
            write(error_unit,'(2A)') "Failed to find name: ", io_exc(i)%shape_name
            write(error_unit,'(A)') "CHECK INPUT"
            STOP
        endif


        !find associated norm entry
        do j=1,size(excitations%norms)
            if(io_exc(i)%norm_name==excitations%norms(j)%name)then
                excitations%exc(i)%int_norm=j
                exit
            endif
        enddo
        if(excitations%exc(i)%int_norm==0)then
            write(error_unit,'(A)') "ERROR, Excitation operator norm name not found in exitation norms"
            write(error_unit,'(2A)') "Failed to find name: ", io_exc(i)%norm_name
            write(error_unit,'(A)') "CHECK INPUT"
            STOP
        endif
    enddo

    !terrible position-get
    i=get_lines('positions.dat')
    allocate(excitations%position(3,i),source=0.0d0)
    call get_position(excitations%position,'positions.dat')

!
!
!    Call set_pos_entry(io,fname,var_name,success)
!    if(.not.success)then
!        write(output_unit,'(/A/)') "No excitations provided in input"
!        return
!    endif
!    read(io,'(a)',iostat=stat)! nothing to read in "excitations" line, one could put the number there
!!   !find out how many entries there are
!
!    nread=0
!    do 
!        read(io,'(a)',iostat=stat) str
!        if (stat < 0)then
!            write(error_unit,'(A)') "end of input file reached reading excitation"
!            exit
!        endif
!        Call io_exc%read_string(str,success)
!        if(success)then
!            nread=nread+1   
!        else
!            exit    !end reached
!        endif
!    enddo
!
!    if(nread>0)then
!        write(output_unit,'(/AI3A/)') "Found ",nread," excitation entries which are read now"
!        !return do beginning of excitation data
!        do i=1,Nread+1
!            backspace(io)
!        enddo
!        excitations%counter=nread
!        allocate(excitations%temporal_param(nread))
!        allocate(excitations%spatial_param (nread))
!        allocate(excitations%name          (nread))
!        do i=1,Nread
!            read(io,'(a)',iostat=stat) str
!            Call io_exc%read_string(str,success)
!            if(.not.success) ERROR STOP "PROGRAMMING MISTAKE, THIS SHOULD ALWAYS WORK"
!            excitations%name(i)=io_exc%excitation_name
!!            excitations%temporal_param(i)%start_value=io_ext%start_value
!!            excitations%temporal_param(i)%end_value  =io_ext%end_value
!!            excitations%temporal_param(i)%start_value=io_ext%start_value
!!            excitations%temporal_param(i)%size_value =io_ext%size_value
!!            excitations%temporal_param(i)%t_start    =io_ext%t_start
!!            excitations%temporal_param(i)%t_end      =io_ext%t_end
!        enddo
!
!    else
!        write(error_unit,'(2/A/A/)') "WARNING, specified excitations, but no excitations are found", "CHECK INPUT"
!    endif
!
!
!!   if(Npair<1)then
!!       write(output_unit,'(/2A/A/)') "Found no entries for ",var_name,' although the keyword is specified'
!!       ERROR STOP "INPUT PROBABLY WRONG"
!!   endif
!!   if(Nnonzero<1)then
!!       write(output_unit,'(/2A/A/)') "Found no nonzero entries for ",var_name,' although the keyword is specified'
!!       ERROR STOP "INPUT PROBABLY WRONG"
!!   endif
!!   write(output_unit,'(/A,I6,2A)') "Found ",Nnonzero," nonzero entries for Hamiltonian ",var_name
!!   !allocate correct size of entries and move IO to beginning of data
!!   allocate(Hpair_tmp(Nnonzero))
!!   do i=1,Npair+1
!!       backspace(io)
!!   enddo
!!   !read in data
!!   ii=1
!!   do i=1,Npair
!!       read(io,*,iostat=stat) attype,dist,val
!!       if(val==0.0d0) cycle
!!       if(attype(2)<attype(1))then
!!           j        =attype(2)
!!           attype(2)=attype(1)
!!           attype(1)=j
!!       endif
!!       Hpair_tmp(ii)%attype=attype
!!       Hpair_tmp(ii)%dist=dist
!!       Hpair_tmp(ii)%val=val
!!       write(output_unit,'(2A,I6,A)') var_name,' entry no.',ii,':'
!!       write(output_unit,'(A,2I6)')    '  atom types:', Hpair_tmp(ii)%attype
!!       write(output_unit,'(A,2I6)')    '  distance  :', Hpair_tmp(ii)%dist
!!       write(output_unit,'(A,E16.8/)') '  energy    :', Hpair_tmp(ii)%val
!!            ii=ii+1
!!        enddo 
end subroutine 

subroutine update_exc(time,lat,dat)
    real(8), intent(in)                     :: time
    type(lattice), intent(inout)            :: lat
    type(excitation_combined),intent(in)    :: dat
    ! internal
    real(8),pointer,contiguous :: opvec(:,:)
    real(8) :: r(3)
    integer :: i,j
    integer :: dim_mode
    integer :: op
    real(8) :: shp(maxval(dim_modes_inner))
  
    do j=1,size(dat%exc)
        op=dat%exc(j)%op
        Call lat%set_order_point_inner(op,opvec)
        dim_mode=dim_modes_inner(op)
        shp(1:dim_mode)=dat%shapes(dat%exc(j)%int_shape)%get_shape(time)
        do i=1,size(dat%position,2)
            r=dat%position(:,i)
            opvec(:,i)= shp(1:dim_mode) * dat%norms(dat%exc(j)%int_norm)%norm(r)
        enddo
    enddo
    nullify(opvec)
end subroutine


!
! routine that chooses if you are doing a EM_cycle, changing temperature or something else
!
subroutine set_excitations(fname,excitation,input_excitations)
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
        if(test)then
            write(*,*) "NAME_VARIABLE ", name_variable
            write(*,*) type_excitations(j)%name
        endif
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
        call read_excitations(io_input,fname,excitations,EM_of_t%temporal_param(j))
        call get_parameter_TPulse(io_input,fname)
    
      case('EMwave')
    
        call read_excitations(io_input,fname,excitations,EM_of_t%temporal_param(j))
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

end subroutine 


subroutine update_EMT_of_r(i,lat)
    implicit none
    integer, intent(in) :: i
    type(lattice), intent(inout) :: lat
    ! internal
    real(8),pointer,contiguous :: opvec(:)
    real(kind=8) :: r(3)
    integer :: j
    integer :: dim_mode
    
    r=position(:,i)
    do j=1,EM_of_t%counter
       Call lat%set_order_point(EM_of_t%temporal_param(j)%op,opvec)
       dim_mode=lat%get_order_dim(EM_of_t%temporal_param(j)%op)
       opvec(1+(i-1)*dim_mode:i*dim_mode)=EM_of_t%spatial_param(j)%field*EM_of_t%spatial_param(j)%norm(r,shape_excitation%center,shape_excitation%cutoff)
    enddo
    STOP  "CONTINUE HERE"
    nullify(opvec)
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
    form=convert('(a,',size_field,'f16.6/)')
    write(6,'(/a,2x,f16.6)') 'real time',time
    write(6,form) 'field value ',EM_of_t%spatial_param(j)%field(:)
  endif
! check data passed incorrectly from dynamics routine
!  if (abs(check(2)).gt.1.0d-8) write(6,'(a,2x,f16.6)') 'Final Temp', check(1)/check(2)/2.0d0/k_B

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
    use m_type_lattice, only: op_name_to_int 
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
            excite%op=op_name_to_int(excite%name)
            if ((excite%name.eq.'Bfield').or.(excite%name.eq.'Efield')) then
                backspace(io)
                read(io,*) dummy,dum_logic,dummy,excite%start_value(1:3),excite%end_value(1:3),excite%t_start,excite%t_end
            else
                backspace(io)
                read(io,*) dummy,dum_logic,dummy,excite%start_value(1),excite%end_value(1),excite%t_start,excite%t_end
            endif
    
       endif
    enddo
end subroutine read_excitations

end module m_excitations
