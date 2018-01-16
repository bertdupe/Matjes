module m_io_utils
interface get_lines
 module procedure get_NB_lines
end interface get_lines

interface get_cols
 module procedure get_NB_columns
end interface get_cols

!interface get_names
! module procedure get_name_simu
!end interface get_names
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The format of the get_parameter subroutine:
! 1. integer - io_unit of the file
! 2. character - name of the file
! 3. character - name of the variable
! 4. variable of the type

interface get_parameter
 module procedure get_2D_vec_real
 module procedure get_1D_vec_real
 module procedure get_1D_vec_int
 module procedure get_1D_vec_bool
 module procedure get_int
 module procedure get_real
 module procedure get_bool
 module procedure get_character
 module procedure get_my_simu
end interface get_parameter

private
public :: get_parameter,get_cols,get_lines

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get the type of simulation that should be done
! inout:  type(type_simu) ALL VARIABLES ARE INITIALIZED TO .FALSE. !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_my_simu(io,my_simu,fname)
use m_derived_types, only : bool_var
use m_simu_parameters, only : type_simu
implicit none
type(bool_var), intent(inout) :: my_simu(:)
integer, intent(in) :: io
character(len=*), intent(in) :: fname
! internal variable
integer :: fin,len_string,nread,i,n_param,ntest,nvar
character(len=30) :: vname
character(len=100) :: str
character(len=10) :: dummy
logical :: success_read

nread=0
ntest=0
nvar=0
success_read=.False.
n_param=size(type_simu)

rewind(io)
do
   read (io,'(a)',iostat=fin) str
   if (fin /= 0) exit
   str= trim(adjustl(str))

   if (len_trim(str)==0) cycle
   if (str(1:1) == '#' ) cycle

   len_string=len_trim(str)
!cccc We start to read the input
   if ( str(1:len_string) == 'simulation') then
      nread=nread+1
      backspace(io)
      read(io,*) dummy, str

      ! find if the variable that was given in input is found in the parameters
      do i=1,n_param
         vname=type_simu(i)%name
         write(*,*) type_simu(i)%name
         ntest=index(vname,str)
          if (ntest.gt.nvar) then
             nvar=ntest
             my_simu(1)%value=.True.
             my_simu(1)%name=type_simu(i)%name
             my_simu(1)%var_name=vname
          endif
      enddo

   endif

enddo

call check_read(nread,my_simu(1)%var_name,fname)

stop

end subroutine get_my_simu

!!! get the names and the values of the type(type_simu)
!subroutine get_name_simu(my_simu,vname)
!use m_derived_types, only : type_simu
!implicit none
!type(type_simu),target, intent(in) :: my_simu
!character(len=*), pointer, intent(out) :: vname(:,:)

!vname(:,1)=my_simu%i_metropolis%name

!end subroutine get_name_simu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! get the name of the variable for derived type
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function get_var_name

!end function get_var_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get a character string (check the string and so on)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_character(io,fname,vname,string)
use m_vector
implicit none
character(len=*), intent(out) :: string
integer, intent(in) :: io
character(len=*), intent(in) :: vname,fname
! internal variable
integer :: fin,len_string,nread
character(len=100) :: str
character(len=10) :: dummy

nread=0
len_string=len(vname)

rewind(io)
do
   read (io,'(a)',iostat=fin) str
   if (fin /= 0) exit
   str= trim(adjustl(str))

   if (len_trim(str)==0) cycle
   if (str(1:1) == '#' ) cycle

!cccc We start to read the input
   if ( str(1:len_string) == vname) then
      nread=nread+1
      backspace(io)
      read(io,*) dummy, string
   endif

enddo

call check_read(nread,vname,fname)

end subroutine get_character

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get a BOOLEAN number parameter (check the string and so on)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_bool(io,fname,vname,tag)
use m_vector
implicit none
logical, intent(out) :: tag
integer, intent(in) :: io
character(len=*), intent(in) :: vname,fname
! internal variable
integer :: fin,len_string,nread
character(len=100) :: str
character(len=10) :: dummy

nread=0
len_string=len(vname)

rewind(io)
do
   read (io,'(a)',iostat=fin) str
   if (fin /= 0) exit
   str= trim(adjustl(str))

   if (len_trim(str)==0) cycle
   if (str(1:1) == '#' ) cycle

!cccc We start to read the input
   if ( str(1:len_string) == vname) then
      nread=nread+1
      backspace(io)
      read(io,*) dummy, tag
   endif

enddo

call check_read(nread,vname,fname)

end subroutine get_bool

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get a REAL number parameter (check the string and so on)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_real(io,fname,vname,number)
use m_vector
implicit none
real(kind=8), intent(out) :: number
integer, intent(in) :: io
character(len=*), intent(in) :: vname,fname
! internal variable
integer :: fin,len_string,nread
character(len=100) :: str
character(len=10) :: dummy

nread=0
len_string=len(vname)

rewind(io)
do
   read (io,'(a)',iostat=fin) str
   if (fin /= 0) exit
   str= trim(adjustl(str))

   if (len_trim(str)==0) cycle
   if (str(1:1) == '#' ) cycle

!cccc We start to read the input
   if ( str(1:len_string) == vname) then
      nread=nread+1
      backspace(io)
      read(io,*) dummy, number
   endif

enddo

call check_read(nread,vname,fname)

end subroutine get_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get a INTEGER number parameter (check the string and so on)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_int(io,fname,vname,number)
use m_vector
implicit none
integer, intent(out) :: number
integer, intent(in) :: io
character(len=*), intent(in) :: vname,fname
! internal variable
integer :: fin,len_string,nread
character(len=100) :: str
character(len=10) :: dummy

nread=0
len_string=len(vname)

rewind(io)
do
   read (io,'(a)',iostat=fin) str
   if (fin /= 0) exit
   str= trim(adjustl(str))

   if (len_trim(str)==0) cycle
   if (str(1:1) == '#' ) cycle

!cccc We start to read the input
   if ( str(1:len_string) == vname) then
      nread=nread+1
      backspace(io)
      read(io,*) dummy, number
   endif

enddo

call check_read(nread,vname,fname)

end subroutine get_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get a BOOLEAN vector parameter (check the string and so on)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_1D_vec_bool(io,fname,vname,N,vec)
use m_vector
implicit none
logical, intent(out) :: vec(:)
integer, intent(in) :: io,N
character(len=*), intent(in) :: vname,fname
! internal variable
integer :: fin,len_string,nread
character(len=100) :: str
character(len=10) :: dummy

nread=0
len_string=len(vname)

rewind(io)
do
   read (io,'(a)',iostat=fin) str
   if (fin /= 0) exit
   str= trim(adjustl(str))

   if (len_trim(str)==0) cycle
   if (str(1:1) == '#' ) cycle

!cccc We start to read the input
   if ( str(1:len_string) == vname) then
      nread=nread+1
      backspace(io)
      read(io,*) dummy, vec(1:N)
   endif

enddo

call check_read(nread,vname,fname)

end subroutine get_1D_vec_bool

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get a REAL vector parameter (check the string and so on)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_1D_vec_real(io,fname,vname,N,vec,vec_norm)
use m_vector
implicit none
real(kind=8), intent(out) :: vec(:)
real(kind=8),optional, intent(in) :: vec_norm
integer, intent(in) :: io,N
character(len=*), intent(in) :: vname,fname
! internal variable
integer :: fin,len_string,nread
character(len=100) :: str
character(len=10) :: dummy
real(kind=8) :: int_norm

nread=0
len_string=len(vname)

rewind(io)
do
   read (io,'(a)',iostat=fin) str
   if (fin /= 0) exit
   str= trim(adjustl(str))

   if (len_trim(str)==0) cycle
   if (str(1:1) == '#' ) cycle

!cccc We start to read the input
   if ( str(1:len_string) == vname) then
      nread=nread+1
      backspace(io)
      read(io,*) dummy, vec(1:N)
   endif

enddo

call check_read(nread,vname,fname)

if (present(vec_norm)) then
! renormalize the variable
   int_norm=norm(vec)
   if (int_norm.gt.1.0d-8) vec=vec/int_norm*vec_norm
endif

end subroutine get_1D_vec_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get a 2D REAL vector parameter (check the string and so on)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_2D_vec_real(io,fname,vname,N,M,vec)
use m_vector
implicit none
real(kind=8), intent(out) :: vec(:,:)
integer, intent(in) :: io,N,M
character(len=*), intent(in) :: vname,fname
! internal variable
integer :: fin,len_string,nread,j
character(len=100) :: str
character(len=10) :: dummy

nread=0
len_string=len(vname)

rewind(io)
do
   read (io,'(a)',iostat=fin) str
   if (fin /= 0) exit
   str= trim(adjustl(str))

   if (len_trim(str)==0) cycle
   if (str(1:1) == '#' ) cycle

!cccc We start to read the input
   if ( str(1:len_string) == vname) then
      nread=nread+1
      backspace(io)
      do j=1,N
         read(io,*) dummy, vec(j,1:M)
      enddo
   endif

enddo

call check_read(nread,vname,fname)

end subroutine get_2D_vec_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get a INTEGER vector parameter (check the string and so on)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_1D_vec_int(io,fname,vname,N,vec)
use m_vector
implicit none
integer, intent(out) :: vec(:)
integer, intent(in) :: io,N
character(len=*), intent(in) :: vname,fname
! internal variable
integer :: fin,len_string,nread
character(len=100) :: str
character(len=10) :: dummy

nread=0
len_string=len(vname)

rewind(io)
do
   read (io,'(a)',iostat=fin) str
   if (fin /= 0) exit
   str= trim(adjustl(str))

   if (len_trim(str)==0) cycle
   if (str(1:1) == '#' ) cycle

!cccc We start to read the input
   if ( str(1:len_string) == vname) then
      nread=nread+1
      backspace(io)
      read(io,*) dummy, vec(1:N)
   endif

enddo

call check_read(nread,vname,fname)

end subroutine get_1D_vec_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get the number of lines in a ASCII file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function get_NB_lines(fname) result(N)
implicit none
character(len=*), intent(in) :: fname
integer :: N
! internal
integer :: nlines,io,fin

nlines=0

open(newunit=io,file=fname,form='formatted',status='old',action='read')
rewind(io)

do
  read (io,*,iostat=fin)
  if (fin/=0) exit
  nlines=nlines+1
enddo
close(io)

N=nlines

end function get_NB_lines

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get the number of columns in a ASCII file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function get_NB_columns(fname) result(N)
implicit none
character(len=*), intent(in) :: fname
integer :: N
! internal
integer :: ncol,io,i
character(len=100) :: str

ncol=0

open(newunit=io,file=fname,form='formatted',status='old',action='read')
rewind(io)

! read the first line and count the columns
read (io,'(a)') str
str= trim(adjustl(str))

close(io)
ncol=0
do i=1,len(str)
if (str(i:i) == '.') ncol=ncol+1
enddo

N=ncol

end function get_NB_columns

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check that the reading went fine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_read(nread,vname,fname)
implicit none
integer, intent(in) :: nread
character(len=*), intent(in) :: vname,fname
! internal variables

select case(nread)
 case(2:)
   write(6,'(3(a,2X))') vname,'was read more than once in file',fname
   stop

 case(0)
   write(6,'(3(a,2X))') vname,'could not be found in file',fname
   stop

 case default
   write(6,'(3(a,2X))') vname,'was read successfully in file',fname

end select

end subroutine check_read

end module m_io_utils
