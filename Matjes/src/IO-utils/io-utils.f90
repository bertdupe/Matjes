module m_io_utils
use m_convert
use m_derived_types

interface get_lines
 module procedure get_NB_lines
end interface get_lines

interface get_cols
 module procedure get_NB_columns
end interface get_cols

interface dump_spinse
 module procedure dump_config_spinse
 module procedure dump_config_spinse_spin
 module procedure dump_config_spinse_vec_point
end interface dump_spinse

interface dump_config
 module procedure dump_config_modes
 module procedure dump_config_matrix_2D_real
 module procedure dump_config_matrix_5D_real
 module procedure dump_config_FFT
 module procedure dump_config_vec_point
 module procedure dump_config_order
end interface dump_config

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
 module procedure get_1D_vec_cmplx
 module procedure get_1D_vec_int
 module procedure get_1D_vec_bool
 module procedure get_int
 module procedure get_real
 module procedure get_bool
 module procedure get_character
 module procedure get_my_simu
 module procedure get_atomic
end interface get_parameter

interface number_nonzero_coeff
  module procedure number_nonzero_coeff_1d
  module procedure number_nonzero_coeff_2d
end interface number_nonzero_coeff


private
public :: get_parameter,get_cols,get_lines,count_variables,get_coeff,dump_config,dump_spinse,number_nonzero_coeff,check_last_char,check_read,max_ind_variable

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! routine that reads and write the Fourrier coefficients layer resolved
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dump_config_FFT(io,kv0,mesh,fftcoef)
use m_vector, only : norm
implicit none
integer, intent(in) :: io
real(kind=8), intent(in) :: kv0(:,:),mesh(:,:,:)
complex(kind=8), intent(in) :: fftcoef(:,:,:)
! internale variables
real(kind=8) :: kk(3)
Integer :: j_lat,i_lat,k
integer :: qnx,qny,Ilat(3),ncoef
character(len=30) :: rw_format

Ilat=shape(fftcoef)
qnx=Ilat(2)
qny=Ilat(3)
ncoef=Ilat(1)

write(rw_format,'( "(", I4, "(2x,f20.15))" )') ncoef+3

do j_lat=1,qny
   do i_lat=1,qnx
      kk=kv0(1,:)*mesh(1,i_lat,j_lat)+kv0(2,:)*mesh(2,i_lat,j_lat)
      write(io,rw_format) kk,(dble(fftcoef(k,i_lat,j_lat)**2+aimag(fftcoef(k,i_lat,j_lat))**2),k=1,ncoef)
    enddo
enddo

end subroutine dump_config_FFT



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! routine that reads and write the local spinse files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dump_config_spinse_spin(io,spin)
use m_constants, only : pi
implicit none
integer, intent(in) :: io
real(kind=8), intent(in) :: spin(:,:)
! internale variables
Integer :: i,shape_spin(2)
real(kind=8) :: widthc,Delta,Bc,Gc,Rc,theta,phi

!     Constants used for the color definition
widthc=5.0d0
Delta =pi*2.0d0/3.0d0
shape_spin=shape(spin)

do i=1,shape_spin(2)

   call get_colors(Rc,Gc,Bc,theta,phi,spin(:,i))

   write(io,'(6(a,f16.8),a)') 'Spin(',theta,',',phi,',',real(i),',',Rc,',',Bc,',',Gc,')'

enddo

end subroutine dump_config_spinse_spin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! routine that reads and write the local spinse files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dump_config_spinse_vec_point(io,spin)
implicit none
integer, intent(in) :: io
type(vec_point), intent(in) :: spin(:)
! internale variables
! internal variables
!     calculating the angles
real(kind=8) :: theta, phi
!     Is the row even (1) or not (0)
Integer :: i,shape_spin
!     colors
real(kind=8) :: Rc,Gc,Bc
!   name of files

shape_spin=size(spin)

do i=1,shape_spin

   call get_colors(Rc,Gc,Bc,theta,phi,spin(i)%w(:))

   write(io,'(6(a,f16.8),a)') 'Spin(',theta,',',phi,',',real(i),',',Rc,',',Bc,',',Gc,')'

enddo

end subroutine dump_config_spinse_vec_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! routine that reads and write the local spinse files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dump_config_spinse(io,my_lattice,position)
use m_derived_types
use m_derived_types, only: lattice
use m_constants, only : pi
implicit none
integer, intent(in) :: io
type(lattice), intent(in) :: my_lattice
real(kind=8), intent(in) :: position(:,:,:,:,:)
! internale variables
Integer :: i_x,i_y,i_z,i_m,N(4)
real(kind=8) :: Rc,Gc,Bc,theta,phi

N=shape(my_lattice%ordpar%l_modes)

do i_m=1,N(4)
   Do i_z=1,N(3)
      Do i_y=1,N(2)
         Do i_x=1,N(1)


        call get_colors(Rc,Gc,Bc,theta,phi,my_lattice%ordpar%l_modes(i_x,i_y,i_z,i_m)%w(:))

         write(io,'(8(a,f16.8),a)') 'Spin(', &
     & theta,',',phi,',',position(1,i_x,i_y,i_z,i_m),',',position(2,i_x,i_y,i_z,i_m),',',position(3,i_x,i_y,i_z,i_m),',', &
     & Rc,',',Bc,',',Gc,')'

         enddo
     enddo
   enddo
enddo

end subroutine dump_config_spinse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! routine that dumps a 5D matrix of real numbers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dump_config_matrix_2D_real(io,matrix)
use m_derived_types
implicit none
integer, intent(in) :: io
real(kind=8), intent(in) :: matrix(:,:)
! internale variables
Integer :: i_x,j,N(2)
character(len=30) :: rw_format

N=shape(matrix)

write(rw_format,'( "(", I4, "f14.8,2x)" )') N(1)

do i_x=1,N(2)

   Write(io,rw_format) (matrix(j,i_x), j=1,N(1))

enddo

end subroutine dump_config_matrix_2D_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! routine that dumps a 5D matrix of real numbers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dump_config_matrix_5D_real(io,matrix)
use m_derived_types
implicit none
integer, intent(in) :: io
real(kind=8), intent(in) :: matrix(:,:,:,:,:)
! internale variables
Integer :: i_x,i_y,i_z,i_m,j_lat,N(5)
character(len=30) :: rw_format

N=shape(matrix)

write(rw_format,'( "(", I4, "f14.8,2x)" )') N(1)

do i_z=1,N(4)
  do i_y=1,N(3)
    do i_x=1,N(2)

    Write(io,rw_format) ((matrix(j_lat,i_x,i_y,i_z,i_m), j_lat=1,N(1)),i_m=1,N(5))

    enddo
  enddo
enddo

end subroutine dump_config_matrix_5D_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! routine that reads and write the local modes configurations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dump_config_modes(io,my_lattice)
use m_derived_types, only : lattice
implicit none
integer, intent(in) :: io
type(lattice), intent(in) :: my_lattice
! internale variables
Integer :: i_x,i_y,i_z,i_m,j_lat,N(4)
character(len=100) :: rw_format

N=shape(my_lattice%ordpar%l_modes)

write(rw_format,'( "(", I4, "f14.8,2x)" )') my_lattice%dim_mode

do i_z=1,N(3)
  do i_y=1,N(2)
    do i_x=1,N(1)

    Write(io,rw_format) ((my_lattice%ordpar%l_modes(i_x,i_y,i_z,i_m)%w(j_lat), j_lat=1,my_lattice%dim_mode),i_m=1,N(4))

    enddo
  enddo
enddo

end subroutine dump_config_modes


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! routine that writes the configuration for order_par
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dump_config_order(io,order)
use m_type_lattice, only : order_par
implicit none
integer, intent(in) :: io
type(order_par), intent(in) :: order
! internale variables
!#Integer :: i_x,i_y,i_z,i_m,j_lat,N(4)
character(len=100) :: rw_format

write(rw_format,'( "(", I4, "f14.8,2x)" )') order%dim_mode
write(io,rw_format) order%all_modes

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! routine that reads and write an array of pointer configurations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dump_config_vec_point(io,modes)
use m_derived_types, only : vec_point
implicit none
integer, intent(in) :: io
type(vec_point), intent(in) :: modes(:)
! internale variables
Integer :: i,N,N_data,j_lat
character(len=20) :: rw_format

N=size(modes)
N_data=size(modes(1)%w)

write(rw_format,'( "(", I4, "f14.8)" )') N_data

do i=1,N

    Write(io,rw_format) (modes(i)%w(j_lat), j_lat=1,N_data)

enddo

end subroutine dump_config_vec_point


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! return the maximal index of parameters of the same type
! in: io tag
! in: variable name (for example J_ D_ or whatever)
! in: name of the file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function max_ind_variable(io,var_name,fname) result(max_ind)
    implicit none
    character(len=*), intent(in) :: var_name,fname
    integer, intent(in) :: io
    integer         ::  max_ind
    ! internal variable
    integer :: length_string,i_var,i_end_str,fin
    character(len=100) :: str
    
    max_ind=0
    length_string=len_trim(var_name)
    
    rewind(io)
    do
        read (io,'(a)',iostat=fin) str
        if (fin /= 0) exit
        str= adjustl(str)

        if (len_trim(str)<length_string) cycle
        if (str(1:1) == '#' ) cycle
        if ( str(1:length_string) == var_name ) then
           i_end_str=scan(str(length_string+1:),' ')
           read(str(length_string+1:length_string+1+i_end_str),*) i_var
           max_ind=max(i_var,max_ind)
        endif
    enddo
    write(6,'(3a,I5)') 'Maximal index of ',var_name,' is:', max_ind
    
end function 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! count the number of parameters of the same type
! in: io tag
! in: varibale name (for example J_ D_ or whatever)
! in: name of the file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer function count_variables(io,var_name,fname)
implicit none
character(len=*), intent(in) :: var_name,fname
integer, intent(in) :: io
! internal variable
integer :: nread,length_string,fin,check,nvariable,n_var_fin
character(len=30) :: str,var_name_local,integer_number
! slope
integer :: i

i=0
nvariable=0
length_string=len_trim(var_name)
n_var_fin=0

rewind(io)
  do
   read (io,'(a)',iostat=fin) str
   if (fin /= 0) exit
   str= trim(adjustl(str))

   if (len_trim(str)==0) cycle
   if (str(1:1) == '#' ) cycle

!cccc We start to read the input
   if ( str(1:length_string) == var_name ) then
      nvariable=nvariable+1
   endif
enddo

! find the largest length of the varibale to read
! this is to differentiate between J_1 and J_10 for example
integer_number=convert(nvariable)
var_name_local=convert(var_name,integer_number)
length_string=len_trim(var_name_local)

do i=1,nvariable
integer_number=convert(i)
var_name_local=convert(var_name,integer_number)
nread=0

  rewind(io)
  do
   read (io,'(a)',iostat=fin) str
   if (fin /= 0) exit
   str= trim(adjustl(str))

   if (len_trim(str)==0) cycle
   if (str(1:1) == '#' ) cycle

!cccc We start to read the input
   if ( str(1:length_string) == var_name_local(1:length_string)) then
      nread=nread+1
      n_var_fin=i
   endif

  enddo

  check=check_read(nread,var_name_local,fname)

enddo

count_variables=n_var_fin

write(6,'(I5,3a)') n_var_fin,' parameters of type ',var_name,' were found'

end function count_variables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get N parameters of the same type
! in: io parameter tag
! in: file name
! in: varibale name (for example J_ D_ or whatever)
! in: matrix of the coefficients
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_coeff(io,fname,var_name,coeff)
implicit none
integer, intent(in) :: io
character(len=*), intent(in) :: var_name,fname
real(kind=8), intent(inout) :: coeff(:)
! internal
integer :: N,i,length_string
character(len=100) :: var_name_local,integer_number

N=size(coeff)
integer_number=convert(N)
var_name_local=convert(var_name,integer_number)
length_string=len_trim(var_name_local)

write(6,'(/)')

do i=1,N
   integer_number=convert(i)
   var_name_local=convert(var_name,integer_number)
   write(*,*) var_name_local
   call get_parameter(io,fname,var_name_local(1:length_string),coeff(i))
enddo

write(6,'(/)')

end subroutine get_coeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get the atomic configuration in the motif
! inout:  type(atom), dimension(:) ALL VARIABLES ARE INITIALIZED TO 0.0 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_atomic(io,fname,var_name,natom,atomic)
use m_basic_types, only : atom
implicit none
type(atom), intent(inout) :: atomic(:)
integer, intent(in) :: io,natom
character(len=*), intent(in) :: var_name,fname
! internal variable
integer :: fin,nread,i,check,length_string,n_read_at
character(len=100) :: str
logical :: success_read

nread=0
n_read_at=0
success_read=.False.
length_string=len_trim(var_name)

rewind(io)
do
   read (io,'(a)',iostat=fin) str
   if (fin /= 0) exit
   str= trim(adjustl(str))

   if (len_trim(str)==0) cycle
   if (str(1:1) == '#' ) cycle

!cccc We start to read the input
   if ( str(1:length_string) == var_name(1:length_string)) then
      nread=nread+1

      ! find if the variable that was given in input is found in the parameters
      do i=1,natom
         n_read_at=n_read_at+1
         read(io,*) atomic(i)%name,atomic(i)%position,atomic(i)%moment
      enddo

   endif

enddo

check=check_read(nread,var_name,fname)

if ((check.eq.0).or.(n_read_at.ne.natom)) stop 'problem reading the unit cell in the unit cell'

end subroutine get_atomic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get the type of simulation that should be done
! inout:  type(type_simu) ALL VARIABLES ARE INITIALIZED TO .FALSE. !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_my_simu(io,fname,my_simu)
use m_derived_types, only : bool_var
use m_simu_parameters, only : type_simu
implicit none
type(bool_var), intent(inout) :: my_simu
integer, intent(in) :: io
character(len=*), intent(in) :: fname
! internal variable
integer :: fin,nread,i,ntest,nvar,n_variable,check
character(len=100) :: vname,vtest
character(len=100) :: str
character(len=100) :: dummy
logical :: success_read

nread=0
ntest=0
nvar=0
success_read=.False.
n_variable=size(type_simu)

rewind(io)
do
   read (io,'(a)',iostat=fin) str
   if (fin /= 0) exit
   str= trim(adjustl(str))

   if (len_trim(str)==0) cycle
   if (str(1:1) == '#' ) cycle

!cccc We start to read the input
   if ( str(1:10) == 'simulation') then
      nread=nread+1
      backspace(io)
      read(io,*) dummy, str
      vtest=trim(str)

      ! find if the variable that was given in input is found in the parameters
      do i=1,n_variable
        vname=type_simu(i)%name
        ntest=index(vname,vtest)

          if (ntest.ne.0) then
             nvar=ntest+nvar
             my_simu%value=.True.
             my_simu%name=type_simu(i)%name
          endif
       if (nvar.eq.1) exit
      enddo

! if the simulation type was not found, write an error message
      if (nvar.eq.0) then
         write(6,'(/,a)') 'The simulation type was not found  '
         write(6,'(2a)') 'The code has read  ', str
         write(6,*) 'possible choices are  ',type_simu
         stop
      endif

   endif

enddo

check=check_read(nread,my_simu%name,fname)

if (check.eq.0) write(6,*) 'default value for variable ', my_simu%name, ' is ', my_simu%value

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
character(len=*), intent(inout) :: string
integer, intent(in) :: io
character(len=*), intent(in) :: vname,fname
! internal variable
integer :: fin,len_string,nread,check
character(len=100) :: str
character(len=100) :: dummy

nread=0
len_string=len(trim(adjustl(vname)))

rewind(io)
do
   read (io,'(a)',iostat=fin) str
   if (fin /= 0) exit
   str= trim(adjustl(str))

   if (len_trim(str)==0) cycle
   if (str(1:1) == '#' ) cycle

!cccc We start to read the input
   if ( str(1:len_string) == trim(adjustl(vname))) then
      nread=nread+1
      backspace(io)
      read(io,*) dummy, string
   endif

enddo

check=check_read(nread,vname,fname)

if (check.eq.0) write(6,*) 'default value for variable ', vname, ' is ', string

end subroutine get_character

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get a BOOLEAN number parameter (check the string and so on)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_bool(io,fname,vname,tag)
use m_vector
implicit none
logical, intent(inout) :: tag
integer, intent(in) :: io
character(len=*), intent(in) :: vname,fname
! internal variable
integer :: fin,len_string,nread,check
character(len=100) :: str
character(len=100) :: dummy

nread=0
len_string=len(trim(adjustl(vname)))

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
      read(io,*) dummy, tag
   endif

enddo

check=check_read(nread,vname,fname)

if (check.eq.0) write(6,*) 'default value for variable ', vname, ' is ', tag

end subroutine get_bool

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get a REAL number parameter (check the string and so on)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_real(io,fname,vname,number)
use m_vector
implicit none
real(kind=8), intent(inout) :: number
integer, intent(in) :: io
character(len=*), intent(in) :: vname,fname
! internal variable
integer :: fin,len_string,nread,check
character(len=100) :: str
character(len=100) :: dummy

nread=0
len_string=len(trim(adjustl(vname)))

rewind(io)
do
   read (io,'(a)',iostat=fin) str
   if (fin /= 0) exit
   str= trim(adjustl(str))

   if (len_trim(str)==0) cycle
   if (str(1:1) == '#' ) cycle

!cccc We start to read the input
   if ( (str(1:len_string) == trim(adjustl(vname))) .and. ( check_last_char(str(len_string+1:len_string+1)) )) then
      nread=nread+1
      backspace(io)
      read(io,*) dummy, number
   endif

enddo

check=check_read(nread,vname,fname)

if (check.eq.0) write(6,*) 'default value for variable ', vname, ' is ', number

end subroutine get_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get a INTEGER number parameter (check the string and so on)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_int(io,fname,vname,number)
use m_vector
implicit none
integer, intent(inout) :: number
integer, intent(in) :: io
character(len=*), intent(in) :: vname,fname
! internal variable
integer :: fin,len_string,nread,check
character(len=100) :: str
character(len=100) :: dummy

nread=0
len_string=len(trim(adjustl(vname)))

rewind(io)
do
   read (io,'(a)',iostat=fin) str
   if (fin /= 0) exit
   str= trim(adjustl(str))

   if (len_trim(str)==0) cycle
   if (str(1:1) == '#' ) cycle

!cccc We start to read the input
   if (( str(1:len_string) == trim(adjustl(vname))) .and. ( check_last_char(str(len_string+1:len_string+1)) )) then
      nread=nread+1
      backspace(io)
      read(io,*) dummy, number
   endif

enddo

check=check_read(nread,vname,fname)

if (check.eq.0) write(6,*) 'default value for variable ', vname, ' is ', number

end subroutine get_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get a BOOLEAN vector parameter (check the string and so on)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_1D_vec_bool(io,fname,vname,N,vec)
use m_vector
implicit none
logical, intent(inout) :: vec(:)
integer, intent(in) :: io,N
character(len=*), intent(in) :: vname,fname
! internal variable
integer :: fin,len_string,nread,check
character(len=100) :: str
character(len=100) :: dummy

nread=0
len_string=len(trim(adjustl(vname)))

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
      read(io,*) dummy, vec(1:N)
   endif

enddo

check=check_read(nread,vname,fname)

if (check.eq.0) write(6,*) 'default value for variable ', vname, ' is ', vec

end subroutine get_1D_vec_bool

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get a REAL vector parameter (check the string and so on)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_1D_vec_real(io,fname,vname,N,vec,vec_norm)
use m_vector
implicit none
real(kind=8), intent(inout) :: vec(:)
real(kind=8),optional, intent(in) :: vec_norm
integer, intent(in) :: io,N
character(len=*), intent(in) :: vname,fname
! internal variable
integer :: fin,len_string,nread,check
character(len=100) :: str
character(len=100) :: dummy
real(kind=8) :: int_norm

nread=0
len_string=len(trim(adjustl(vname)))

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
      read(io,*) dummy, vec(1:N)
   endif

enddo

check=check_read(nread,vname,fname)

if (present(vec_norm)) then
! renormalize the variable
   int_norm=norm(vec)
   if (int_norm.gt.1.0d-8) vec=vec/int_norm*vec_norm
endif

if (check.eq.0) write(6,*) 'default value for variable ', vname, ' is ', vec

end subroutine get_1D_vec_real



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get a Complex vector parameter (check the string and so on)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_1D_vec_cmplx(io,fname,vname,N,vec)
use m_vector
implicit none
complex(kind=8), intent(inout) :: vec(:)
integer, intent(in) :: io,N
character(len=*), intent(in) :: vname,fname
! internal variable
integer :: fin,len_string,nread,check
character(len=100) :: str
character(len=100) :: dummy

nread=0
len_string=len(trim(adjustl(vname)))

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
      read(io,*) dummy, vec(1:N)
   endif

enddo

check=check_read(nread,vname,fname)

if (check.eq.0) write(6,*) 'default value for variable ', vname, ' is ', vec

end subroutine get_1D_vec_cmplx


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get a 2D REAL vector parameter (check the string and so on)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_2D_vec_real(io,fname,vname,N,M,vec)
use m_vector
implicit none
real(kind=8), intent(inout) :: vec(:,:)
integer, intent(in) :: io,N,M
character(len=*), intent(in) :: vname,fname
! internal variable
integer :: fin,len_string,nread,j,check
character(len=100) :: str

nread=0
len_string=len(trim(adjustl(vname)))

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
      do j=1,N
         read(io,*) vec(j,1:M)
      enddo
   endif

enddo

check=check_read(nread,vname,fname)

if (check.eq.0) write(6,*) 'default value for variable ', vname, ' is ', vec

end subroutine get_2D_vec_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get a INTEGER vector parameter (check the string and so on)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_1D_vec_int(io,fname,vname,N,vec)
use m_vector
implicit none
integer, intent(inout) :: vec(:)
integer, intent(in) :: io,N
character(len=*), intent(in) :: vname,fname
! internal variable
integer :: fin,len_string,nread,check
character(len=100) :: str
character(len=100) :: dummy

nread=0
len_string=len(trim(adjustl(vname)))

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
      read(io,*) dummy, vec(1:N)
   endif

enddo

check=check_read(nread,vname,fname)

if (check.eq.0) write(6,*) 'default value for variable ', vname, ' is ', vec

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
character(len=400) :: str
logical :: test

ncol=0
N=0
inquire(FILE=fname,EXIST=test)
if (.not.test) then
  write(6,'(a,2x,a,2x,a)') 'file',fname,'not found'
  return
endif

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
! check the last character of a string (should be ' ' or '=')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
logical function check_last_char(string)
implicit none
character(len=1), intent(in) :: string

check_last_char=.false.

if ( (string.eq.' ') .or. (string.eq.'=') ) check_last_char=.true.

end function check_last_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check that the reading went fine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function check_read(nread,vname,fname)
implicit none
integer, intent(in) :: nread
character(len=*), intent(in) :: vname,fname
! internal variables
check_read=-10

select case(nread)
 case(2:)
   write(6,'(3(a,2X))') vname,'was read more than once in file ',fname
   stop

 case(0)
   write(6,'(3(a,2X))') vname,'could not be found in file ',fname
   check_read=0

 case default
   write(6,'(3(a,2X))') vname,'was read successfully in file ',fname
   check_read=1

end select

if (check_read.eq.-10) then
   write(6,'(2(a,2X))') 'problem reading the variable ', vname
   stop
endif

end function check_read

! ===============================================================
! function that gets the RGB colors
! ===============================================================
subroutine get_colors(Rc,Gc,Bc,theta,phi,mode)
use m_constants,only : pi
implicit none
real(kind=8),intent(inout) :: Rc,Gc,Bc,theta,phi
real(kind=8),intent(in) :: mode(3)
! internal
real(kind=8) :: widthc,Delta,phi_color

widthc=5.0d0
Delta =pi*2.0d0/3.0d0

!       Yes for these formulars it is helpfull to make a picture.
!       The initial object is a cone with its top at
!       coordinates (0,0,1). First turn it around the y-achse into
!       the x-z-plane around angly, then turn it around the z-achse
!       into the right position.
!       Then translate it to the right r_x,r_y position.

if (abs(mode(3)).lt.1.0d0) then
  theta=acos(mode(3))*180.0d0/pi
else
  theta=90.0d0-dsign(90.0d0,mode(3))
endif

phi=atan2(mode(2),mode(1))

phi=phi*180.0d0/pi

!       Calcualting the color as a function of the angle in or
!       out of the plane
phi_color=pi*theta/300.0d0*2.0d0
Rc = widthc*(cos(phi_color+0*Delta))
if (Rc.lt.0.000001d0)  Rc=0.0d0
Gc = widthc*(cos(phi_color+1*Delta))
if (Gc.lt.0.000001d0)  Gc=0.0d0
Bc = widthc*(cos(phi_color+2*Delta))
if (Bc.lt.0.000001d0)  Bc=0.0d0

end subroutine get_colors
! ===============================================================



!
! function that determines how many coefficients are none 0 in the exchange or the DMI or whatever
!
integer function number_nonzero_coeff_2d(coeff,name)
implicit none
real(kind=8), intent(in) :: coeff(:,:)
character(len=*), intent(in) :: name
! internal
integer :: i,N(2)
real(kind=8) :: norm

N=shape(coeff(:,:))
number_nonzero_coeff_2d=0

! loop over the number of coefficients
do i=1,N(2)
   norm=sqrt(coeff(1,i)**2+coeff(2,i)**2+coeff(3,i)**2)
   if (norm.gt.1.0d-8) number_nonzero_coeff_2d=i
enddo

if (number_nonzero_coeff_2d.eq.0) then
   write(6,'(2a)') 'no non-zero coefficients found for ',name
endif

end function number_nonzero_coeff_2d

integer function number_nonzero_coeff_1d(coeff,name)
implicit none
real(kind=8), intent(in) :: coeff(:)
character(len=*), intent(in) :: name
! internal
integer :: i,N

N=size(coeff)
number_nonzero_coeff_1d=0

! loop over the number of coefficients
do i=1,N
   if (abs(coeff(i)).gt.1.0d-8) number_nonzero_coeff_1d=i
enddo

if (number_nonzero_coeff_1d.eq.0) then
   write(6,'(2a)') 'no non-zero coefficients found for ',name
endif

end function number_nonzero_coeff_1d

end module m_io_utils
