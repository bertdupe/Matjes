module m_io_utils
use m_convert
use m_derived_types
use m_io_read_util
use, intrinsic :: iso_fortran_env, only : output_unit,error_unit
implicit none

interface get_lines
 module procedure get_NB_lines
end interface get_lines

interface get_cols
 module procedure get_NB_columns
end interface get_cols

interface dump_spinse
 module procedure dump_config_spinse
 module procedure dump_config_spinse_spin
end interface dump_spinse

interface dump_config
 module procedure dump_config_modes
 module procedure dump_config_matrix_N_1D_real
 module procedure dump_config_matrix_2D_real
 module procedure dump_config_matrix_5D_real
 module procedure dump_config_FFT
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
 module procedure get_vec1D_real
 module procedure get_vec1D_int
 module procedure get_int
 module procedure get_int8
 module procedure get_real
 module procedure get_bool
 module procedure get_character
 module procedure get_my_simu
 module procedure get_atomic
 module procedure get_atomtype
 module procedure get_cell
 module procedure get_H_pair
 module procedure get_H_pair_tensor
 module procedure get_H_triple
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
Integer :: i,shape_spin(2),j,counter
real(kind=8) :: widthc,Delta,Bc,Gc,Rc,theta,phi

!     Constants used for the color definition
widthc=5.0d0
Delta =pi*2.0d0/3.0d0
shape_spin=shape(spin)
counter=0

do i=1,shape_spin(2)
  do j=1,shape_spin(1),3

   counter=counter+1

   call get_colors(Rc,Gc,Bc,theta,phi,spin(j:j+2,i))

   write(io,'(6(a,f16.8),a)') 'Spin(',theta,',',phi,',',real(counter),',',Rc,',',Bc,',',Gc,')'
  enddo
enddo

end subroutine dump_config_spinse_spin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! routine that reads and write the local spinse files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dump_config_spinse(io,lat,position)
!    use m_derived_types
    use m_derived_types, only: lattice
    use m_constants, only : pi
    implicit none
    integer, intent(in) :: io
    type(lattice), intent(in) :: lat
    real(kind=8), intent(in) :: position(:,:,:,:,:)
    ! internale variables
    Integer :: i_x,i_y,i_z,i_m
    real(kind=8) :: Rc,Gc,Bc,theta,phi


    do i_m=1,lat%M%dim_mode/3
       Do i_z=1,lat%dim_lat(3)
          Do i_y=1,lat%dim_lat(2)
             Do i_x=1,lat%dim_lat(1)
    
            call get_colors(Rc,Gc,Bc,theta,phi,lat%M%modes((i_m-1)*3+1:i_m*3,i_x,i_y,i_z))
    
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
! routine that dumps a 1D matrix of real numbers with N entries per line
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dump_config_matrix_N_1D_real(io,N,matrix)
    implicit none
    integer, intent(in) :: io
    integer, intent(in) :: N    !
    real(8), intent(in) :: matrix(:)
    ! internale variables
    character(len=30) :: rw_format
    
    write(rw_format,'( "(", I4, "E24.16,2x)" )') N
    write(io,rw_format) matrix
end subroutine 


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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! routine that reads and write the local modes configurations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dump_config_modes(io,lat)
!this probably is supposed to write out all modes, but right now I updated it to
!just write out the M mode since I am not sure if this is still used
use m_derived_types, only : lattice
implicit none
integer, intent(in) :: io
type(lattice), intent(in) :: lat
! internale variables
Integer :: i_x,i_y,i_z,i_m,j_lat,N(4)
character(len=100) :: rw_format


N(1:3)=lat%dim_lat

write(rw_format,'( "(", I4, "E24.16,2x)" )') lat%M%dim_mode ! have all the modes of the unit cell on one line

do i_z=1,N(3)
  do i_y=1,N(2)
    do i_x=1,N(1)

    Write(io,rw_format) ((lat%M%modes(i_m,i_x,i_y,i_z), j_lat=1,lat%M%dim_mode),i_m=1,N(4))

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

write(rw_format,'( "(", I4, "f14.8,2x)" )') order%dim_mode ! this is to have the all the order parameter of a unit cell on a line
!write(rw_format,'( "(", I4, "f14.8,2x)" )') 3   ! one order parameter per line.

write(io,rw_format) order%all_modes

end subroutine


function max_ind_variable(io,var_name,fname) result(max_ind)
    implicit none
    character(len=*), intent(in) :: var_name,fname
    integer, intent(in) :: io
    integer         ::  max_ind
    ! internal variable
    integer :: length_string,i_var,i_end_str,fin,stat
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
           read(str(length_string+1:length_string+1+i_end_str),*,iostat=stat) i_var
           if(stat==0) max_ind=max(i_var,max_ind)
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
! in: stride_in, optional striding parameter which allows read-in of array values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_coeff(io,fname,var_name,coeff,stride_in)
implicit none
integer, intent(in) :: io
character(len=*), intent(in) :: var_name,fname
real(kind=8), intent(inout) :: coeff(:)
integer,intent(in),optional :: stride_in
! internal
integer :: N,i,length_string
character(len=100) :: var_name_local,integer_number
integer :: stride

stride=1
if(present(stride_in)) stride=stride_in

N=size(coeff)/stride
integer_number=convert(N)
var_name_local=convert(var_name,integer_number)
length_string=len_trim(var_name_local)

write(6,'(/)')

do i=1,N
   integer_number=convert(i)
   var_name_local=convert(var_name,integer_number)
   length_string=len_trim(var_name_local)
   call get_parameter(io,fname,var_name_local(1:length_string),stride,coeff((i-1)*stride+1:i*stride))
enddo

write(6,'(/)')

end subroutine get_coeff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get cell based on atomtype
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine get_H_triple(io,fname,var_name,Htriples,success)
    use m_input_H_types, only: Hr_triple,reduce_Hr_triple,Hr_triple_single
    !always uses the same variable name, hence kept as parameter

    integer, intent(in)                         :: io
    character(len=*), intent(in)                :: fname,var_name
    type(Hr_triple), intent(out), allocatable   :: Htriples(:)
    logical,intent(out)                         :: success
    ! internal variable
    type(Hr_triple_single),allocatable   :: Htriple_tmp(:)
    type(Hr_triple), allocatable         :: Htriples_tmp(:)
    integer :: attype(3),dist
    real(8) :: val

    integer :: Ntriple,Nnonzero
    integer :: nread,i,ii,j
    integer :: stat
    character(len=100) :: str

    nread=0
    Call set_pos_entry(io,fname,var_name,success)
    if(success)then
        success=.false.
        read(io,'(a)',iostat=stat) str
        Ntriple=0; Nnonzero=0
        nread=nread+1
        do 
            read(io,'(a)',iostat=stat) str
            if (stat /= 0) STOP "UNEXPECTED END OF INPUT FILE"
            read(str,*,iostat=stat) attype,dist,val
            if (stat /= 0) exit
            Ntriple=Ntriple+1
            if(val/=0.0d0) Nnonzero=Nnonzero+1
        enddo
        if(Ntriple<1)then
            write(error_unit,'(/2A/A/)') "Found no entries for ",var_name,' although the keyword is specified'
#ifndef CPP_SCRIPT            
            ERROR STOP "INPUT PROBABLY WRONG (disable with CPP_SCRIPT preprocessor flag)"
#endif
            return
        endif
        if(Nnonzero<1)then
            write(error_unit,'(/2A/A/)') "Found no nonzero entries for ",var_name,' although the keyword is specified'
#ifndef CPP_SCRIPT            
            ERROR STOP "INPUT PROBABLY WRONG (disable with CPP_SCRIPT preprocessor flag)"
#endif
            return
        endif
        write(output_unit,'(/A,I6,2A)') "Found ",Nnonzero," nonzero entries for Hamiltonian ",var_name
        success=.true.
        allocate(Htriple_tmp(Nnonzero))
        do i=1,Ntriple+1
            backspace(io)
        enddo
        ii=1
        do i=1,Ntriple
            read(io,*,iostat=stat) attype,dist,val
            if(val==0.0d0) cycle
            if(attype(2)<attype(1))then
                j        =attype(2)
                attype(2)=attype(1)
                attype(1)=j
            endif
            Htriple_tmp(ii)%attype=attype
            Htriple_tmp(ii)%dist=dist
            Htriple_tmp(ii)%val=val
            write(output_unit,'(2A,I6,A)') var_name,' entry no.',ii,':'
            write(output_unit,'(A,3I6)')    '  atom types:', Htriple_tmp(ii)%attype
            write(output_unit,'(A,2I6)')    '  distance  :', Htriple_tmp(ii)%dist
            write(output_unit,'(A,E16.8/)') '  energy    :', Htriple_tmp(ii)%val
            ii=ii+1
        enddo 

        !combines single entries into arrays with same atom types
        Call reduce_Hr_triple(Htriple_tmp,Htriples) 
        !check if any entry appears more than once
        do i=1,size(Htriples)
            do j=2,size(Htriples(i)%dist)
                if(any(Htriples(i)%dist(j)==Htriples(i)%dist(:j-1)))then
                    write(output_unit,*) "ERROR, found the same distance twice for for Hamiltonian ",var_name
                    write(output_unit,*) "       at atom indices  :",Htriples(i)%attype
                    write(output_unit,*) "       with the distance:",Htriples(i)%dist(j)
                    STOP "THIS IS MOST PROBABLY AN INPUT MISTAKE"
                endif
            enddo
        enddo

        !symmetrize different type Hamiltonians  (i.e. all attype=[1 2] interactions are dublicated with [2 1]
        if(any(Htriples%attype(1)/=Htriples%attype(2)))then
            Call move_alloc(Htriples,Htriples_tmp)
            allocate(Htriples(size(Htriples_tmp)+count(Htriples_tmp%attype(1)/=Htriples_tmp%attype(2))))
            ii=0
            do i=1,size(Htriples_tmp)
                ii=ii+1
                Htriples(ii)=Htriples_tmp(i)
                if(Htriples_tmp(i)%attype(1)/=Htriples_tmp(i)%attype(2))then
                    ii=ii+1
                    Htriples(ii)=Htriples_tmp(i)
                    Htriples(ii)%attype(1:2)=[Htriples_tmp(i)%attype(2),Htriples_tmp(i)%attype(1)]
                endif
            enddo
            deallocate(Htriples_tmp)
        endif

        success=.true.
    endif

    Call check_further_entry(io,fname,var_name)
end subroutine 



subroutine get_H_pair(io,fname,var_name,Hpairs,success)
    use m_input_H_types, only: Hr_pair,reduce_Hr_pair,Hr_pair_single
    !always uses the same variable name, hence kept as parameter

    integer, intent(in)                         :: io
    character(len=*), intent(in)                :: fname,var_name
    type(Hr_pair), intent(out), allocatable     :: Hpairs(:)
    logical,intent(out)                         :: success
    ! internal variable
    type(Hr_pair_single),allocatable            :: Hpair_tmp(:)
    type(Hr_pair), allocatable                  :: Hpairs_tmp(:)
    integer :: attype(2),dist
    real(8) :: val

    integer :: Npair,Nnonzero
    integer :: nread,i,ii,j
    integer :: stat
    character(len=100) :: str
    character(len=100) :: comment_char, dummy

    nread=0
    Call set_pos_entry(io,fname,var_name,success)
    read(io,'(a)',iostat=stat) str
    if(success)then
        !find out how many entries there are
        success=.false.
        Npair=0; Nnonzero=0
        nread=nread+1
        do 
            read(io,'(a)',iostat=stat) str
            if (stat /= 0) STOP "UNEXPECTED END OF INPUT FILE"
            read(str,*,iostat=stat) comment_char, dummy
            if (comment_char(1:1) .eq. '#') cycle
            read(str,*,iostat=stat) attype,dist,val
            if (stat /= 0) exit
            Npair=Npair+1
            if(val/=0.0d0) Nnonzero=Nnonzero+1
        enddo
        if(Npair<1)then
            write(error_unit,'(/2A/A/)') "Found no entries for ",var_name,' although the keyword is specified'
#ifndef CPP_SCRIPT            
            ERROR STOP "INPUT PROBABLY WRONG (disable with CPP_SCRIPT preprocessor flag)"
#endif
            return
        endif
        if(Nnonzero<1)then
            write(error_unit,'(/2A/A/)') "WARNING, Found no nonzero entries for: ",var_name,' although the keyword is specified'
#ifndef CPP_SCRIPT            
            ERROR STOP "INPUT PROBABLY WRONG (disable with CPP_SCRIPT preprocessor flag)"
#endif
            return
        endif
        write(output_unit,'(/A,I6,2A)') "Found ",Nnonzero," nonzero entries for Hamiltonian ",var_name
        !allocate correct size of entries and move IO to beginning of data
        allocate(Hpair_tmp(Nnonzero))
        do i=1,Npair+1
            backspace(io)
        enddo
        !read in data
        ii=1
        do i=1,Npair
            read(io,*,iostat=stat) attype,dist,val
            if(val==0.0d0) cycle
            if(attype(2)<attype(1))then
                j        =attype(2)
                attype(2)=attype(1)
                attype(1)=j
            endif
            Hpair_tmp(ii)%attype=attype
            Hpair_tmp(ii)%dist=dist
            Hpair_tmp(ii)%val=val
            write(output_unit,'(2A,I6,A)') var_name,' entry no.',ii,':'
            write(output_unit,'(A,2I6)')    '  atom types:', Hpair_tmp(ii)%attype
            write(output_unit,'(A,2I6)')    '  distance  :', Hpair_tmp(ii)%dist
            write(output_unit,'(A,E16.8/)') '  energy    :', Hpair_tmp(ii)%val
            ii=ii+1
        enddo 

        !combines single entries into arrays with same atom types
        Call reduce_Hr_pair(Hpair_tmp,Hpairs) 
        !check if any entry appears more than once
        do i=1,size(Hpairs)
            do j=2,size(Hpairs(i)%dist)
                if(any(Hpairs(i)%dist(j)==Hpairs(i)%dist(:j-1)))then
                    write(output_unit,*) "ERROR, found the same distance twice for Hamiltonian ",var_name
                    write(output_unit,*) "       at atom indices  :",Hpairs(i)%attype
                    write(output_unit,*) "       with the distance:",Hpairs(i)%dist(j)
                    STOP "THIS IS MOST PROBABLY AN INPUT MISTAKE"
                endif
            enddo
        enddo

        !symmetrize different type Hamiltonians  (i.e. all attype=[1 2] interactions are dublicated with [2 1]
        if(any(Hpairs%attype(1)/=Hpairs%attype(2)))then
            Call move_alloc(Hpairs,Hpairs_tmp)
            allocate(Hpairs(size(Hpairs_tmp)+count(Hpairs_tmp%attype(1)/=Hpairs_tmp%attype(2))))
            ii=0
            do i=1,size(Hpairs_tmp)
                ii=ii+1
                Hpairs(ii)=Hpairs_tmp(i)
                if(Hpairs_tmp(i)%attype(1)/=Hpairs_tmp(i)%attype(2))then
                    ii=ii+1
                    Hpairs(ii)=Hpairs_tmp(i)
                    Hpairs(ii)%attype=[Hpairs_tmp(i)%attype(2),Hpairs_tmp(i)%attype(1)]
                endif
            enddo
            deallocate(Hpairs_tmp)
        endif

        success=.true.
    endif

    Call check_further_entry(io,fname,var_name)
end subroutine 

subroutine get_H_pair_tensor(io,fname,var_name,Hpairs_tensor,success)
    use m_input_H_types, only: reduce_Hr_tensor_pair,Hr_pair_tensor,Hr_pair_single_tensor
    !always uses the same variable name, hence kept as parameter

    integer, intent(in)                            :: io
    character(len=*), intent(in)                   :: fname,var_name
    type(Hr_pair_tensor), intent(out), allocatable :: Hpairs_tensor(:)
    logical,intent(out)                            :: success
    ! internal variable
    type(Hr_pair_single_tensor),allocatable        :: Hpair_tmp(:)
    type(Hr_pair_tensor), allocatable              :: Hpairs_tensor_tmp(:)
    integer :: attype(2),dist
    real(8) :: val(9),bound(3)

    integer :: Npair,Nnonzero
    integer :: nread,i,ii,j
    integer :: stat
    character(len=100) :: str

    nread=0
    val=0.0d0
    bound=0.0d0
    Call set_pos_entry(io,fname,var_name,success)
    read(io,'(a)',iostat=stat) str
    if(success)then
        !find out how many entries there are
        success=.false.
        Npair=0; Nnonzero=0
        nread=nread+1
        do
            read(io,'(a)',iostat=stat) str
            if (stat /= 0) STOP "UNEXPECTED END OF INPUT FILE"
            read(str,*,iostat=stat) attype,dist,val
            if (stat /= 0) exit
            Npair=Npair+1
            do i=1,9
               if(val(i)/=0.0d0) then
                   Nnonzero=Nnonzero+1
                   exit
               endif
            enddo
        enddo
        if(Npair<1)then
            write(error_unit,'(/2A/A/)') "Found no entries for ",var_name,' although the keyword is specified'
#ifndef CPP_SCRIPT
            ERROR STOP "INPUT PROBABLY WRONG (disable with CPP_SCRIPT preprocessor flag)"
#endif
            return
        endif
        if(Nnonzero<1)then
            write(error_unit,'(/2A/A/)') "WARNING, Found no nonzero entries for: ",var_name,' although the keyword is specified'
#ifndef CPP_SCRIPT
            ERROR STOP "INPUT PROBABLY WRONG (disable with CPP_SCRIPT preprocessor flag)"
#endif
            return
        endif
        write(output_unit,'(/A,I6,2A)') "Found ",Nnonzero," nonzero entries for Hamiltonian ",var_name
        !allocate correct size of entries and move IO to beginning of data
        allocate(Hpair_tmp(Nnonzero))
        do i=1,Npair+1
            backspace(io)
        enddo
        !read in data
        ii=1
        do i=1,Npair
            val=0.0d0
            read(io,*,iostat=stat) attype,dist,val,bound
            if(all(val==0.0d0)) cycle
            if(attype(2)<attype(1))then
                j        =attype(2)
                attype(2)=attype(1)
                attype(1)=j
            endif
            Hpair_tmp(ii)%attype=attype
            Hpair_tmp(ii)%dist=dist
            Hpair_tmp(ii)%val=val
            Hpair_tmp(ii)%bound=bound
            write(output_unit,'(2A,I6,A)') var_name,' entry no.',ii,':'
            write(output_unit,'(A,2I6)')    '  atom types:', Hpair_tmp(ii)%attype
            write(output_unit,'(A,2I6)')    '  distance  :', Hpair_tmp(ii)%dist
            write(output_unit,'(A,9E16.8)') '  energy    :', Hpair_tmp(ii)%val
            write(output_unit,'(A,3E16.8/)') ' along the bound :', Hpair_tmp(ii)%bound
            ii=ii+1
        enddo

        !combines single entries into arrays with same atom types
        Call reduce_Hr_tensor_pair(Hpair_tmp,Hpairs_tensor)

        !check if any entry appears more than once

        do i=1,size(Hpairs_tensor)
            do j=2,size(Hpairs_tensor(i)%dist)
                if(any(Hpairs_tensor(i)%dist(j)==Hpairs_tensor(i)%dist(:j-1)))then
                    write(output_unit,*) "ERROR, found the same distance twice for Hamiltonian ",var_name
                    write(output_unit,*) "       at atom indices  :",Hpairs_tensor(i)%attype
                    write(output_unit,*) "       with the distance:",Hpairs_tensor(i)%dist(j)
                    STOP "THIS IS MOST PROBABLY AN INPUT MISTAKE"
                endif
            enddo
        enddo

        !symmetrize different type Hamiltonians  (i.e. all attype=[1 2] interactions are duplicated with [2 1]
        if(any(Hpairs_tensor%attype(1)/=Hpairs_tensor%attype(2)))then
            Call move_alloc(Hpairs_tensor,Hpairs_tensor_tmp)
            allocate(Hpairs_tensor(size(Hpairs_tensor_tmp)+count(Hpairs_tensor_tmp%attype(1)/=Hpairs_tensor_tmp%attype(2))))
            ii=0
            do i=1,size(Hpairs_tensor_tmp)
                ii=ii+1
                Hpairs_tensor(ii)=Hpairs_tensor_tmp(i)
                if(Hpairs_tensor_tmp(i)%attype(1)/=Hpairs_tensor_tmp(i)%attype(2))then
                    ii=ii+1
                    Hpairs_tensor(ii)=Hpairs_tensor_tmp(i)
                    Hpairs_tensor(ii)%attype=[Hpairs_tensor_tmp(i)%attype(2),Hpairs_tensor_tmp(i)%attype(1)]
                endif
            enddo
            deallocate(Hpairs_tensor_tmp)
        endif

        success=.true.
    endif

    Call check_further_entry(io,fname,var_name)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get cell based on atomtype
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_cell(io,fname,attype,cell)
    !always uses the same variable name, hence kept as parameter

    integer, intent(in)             :: io
    character(len=*), intent(in)    :: fname
    type(atomtype), intent(in)      :: attype(:)
    type(t_cell),intent(out)        :: cell
    ! internal variable
    integer     ::  N_atoms
    character(len=*),parameter  ::  var_name='atoms'
    character(len=30)       ::  atname
    integer     :: id
    integer :: i,j,length_string
    integer :: stat
    character(len=100) :: str
    logical :: used(size(attype))
    integer :: div

    cell%n_attype=size(attype)
    Call set_pos_entry(io,fname,var_name)
    length_string=len_trim(var_name)
    read (io,'(a)',iostat=stat) str
    read(str(length_string+1:),*,iostat=stat) N_atoms
    if(stat/=0) STOP "atoms keyword found in input, but the number of atoms is not specified there"
    write(output_unit,'(/A,I6,A)') "found atomcell keyword, start reading the ",N_atoms," atoms"
    allocate(cell%atomic(N_atoms))
    do i=1,N_atoms
        read(io,'(a)',iostat=stat) str
        if (stat /= 0) STOP "UNEXPECTED END OF INPUT FILE"
        read(str,*,iostat=stat) atname,cell%atomic(i)%position, div
        if (stat == 0)then
            cell%atomic(i)%position=cell%atomic(i)%position/real(div,8)
        else
            read(str,*,iostat=stat) atname,cell%atomic(i)%position
            if (stat /= 0) STOP "FAILED TO READ ATOMIC id/name and position"
        endif
        read(atname,*,iostat=stat) id
        if(stat/=0)then
            id=0
            do j=1,size(attype)
                if(trim(attype(j)%name)==trim(atname))then
                    id=j
                    exit
                endif
            enddo
            if(id==0) STOP "Did not find name of atom in atom types"
        endif
        if(id<1.or.id>size(attype)) STOP "ATOM TYPE ID IS NONSENSICAL"
        cell%atomic(i)%type_id=id
        Call cell%atomic(i)%set_attype(attype(id))
        write(output_unit,*) "successfully read atom number:",i," of ",N_atoms
    enddo

    Call check_further_entry(io,fname,var_name)

    !check all types are used
    used=[(any(cell%atomic(:)%type_id==i),i=1,size(attype))]
    do j=1,size(used)
        if(.not.used(j)) write(error_unit,'(3A)') 'ATOM TYPE "',trim(attype(j)%name), '" is not used'
    enddo
    if(.not.all(used)) STOP "SOME ATOM TYPES ARE NOT USED, CHECK CELL SETUP"
end subroutine 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get the different atomic types used in the unit-cell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_atomtype(io,fname,attype)
    !always uses the same variable name, hence kept as parameter
    use m_basic_types, only : atom

    implicit none
    type(atomtype), intent(inout),allocatable :: attype(:)
    integer, intent(in) :: io
    character(len=*), intent(in) :: fname
    ! internal variable
    integer     ::  N_attype,j
    character(len=*),parameter  ::  var_name='atomtypes'
    integer :: i,length_string
    integer :: stat
    character(len=100) :: str

    Call set_pos_entry(io,fname,var_name)
    
    length_string=len_trim(var_name)
    read (io,'(a)',iostat=stat) str
    read(str(length_string+1:),*,iostat=stat) N_attype
    if(stat/=0) STOP "atomtypes keyword found in input, but the number of atom types is not specified there"
    if(allocated(attype)) STOP "PROGRAMMING MISTAKE? attype should not be allocated at start of get_atomic"
    allocate(attype(N_attype))
    do i=1,N_attype
        write(output_unit,*) "Trying to read atom type ",i," of ",N_attype
        read(io,'(a)',iostat=stat) str
        if (stat /= 0) STOP "UNEXPECTED END OF INPUT FILE"
        read(str,*,iostat=stat) attype(i)%name,attype(i)%moment,attype(i)%charge,attype(i)%mass,attype(i)%use_ph,attype(i)%orbitals
        if (stat == 0)then
            write(output_unit,*) attype(i)
            cycle
        endif
        read(str,*,iostat=stat) attype(i)%name,attype(i)%moment,attype(i)%charge,attype(i)%mass,attype(i)%use_ph
        if (stat == 0)then
            write(output_unit,*) attype(i)
            cycle
        endif
        STOP "FAILED TO READ PARAMETERS IN ATOMTYPE, ARE ALL RELEVANT PARAMETERS SPECIFIED AND ARE ALL ATOM TYPES SET?"
    enddo

    Call check_further_entry(io,fname,var_name)

    do j=1,N_attype
        do i=1,j-1
            if(trim(attype(i)%name)==trim(attype(j)%name))then
                write(error_unit,*) "Error, found 2 atoms types with the same name: atom type",i," and ", j
                write(error_unit,*) "atom type ",i
                write(error_unit,*) attype(i)
                write(error_unit,*) 
                write(error_unit,*) "atom type ",j
                write(error_unit,*) attype(j)
                STOP
            endif
        enddo
    enddo
end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get the atomic configuration in the motif
! inout:  type(atom), dimension(:) ALL VARIABLES ARE INITIALIZED TO 0.0 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_atomic(io,fname,var_name,natom,atomic)
!old and obsolete way to read the atomic input in
use m_basic_types, only : atom
implicit none
type(atom), intent(inout) :: atomic(:)
integer, intent(in) :: io,natom
character(len=*), intent(in) :: var_name,fname
! internal variable
integer :: fin,nread,i,check,length_string,n_read_at
integer :: stat
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
         atomic(i)%type_id=n_read_at    !with old input all atoms are of different type
         read(io,'(a)',iostat=fin) str
         read(str,*,iostat=stat) atomic(i)%name,atomic(i)%position,atomic(i)%moment,atomic(i)%charge
         if(stat==0)then
            write(*,*) "Read magnetic moment and charge from atom",i
         else
            read(str,*,iostat=stat) atomic(i)%name,atomic(i)%position,atomic(i)%moment
            if(stat==0)then
                write(*,*) "Read magnetic moment from atom",i
            else
                read(str,*,iostat=stat) atomic(i)%name,atomic(i)%position
                if(stat/=0)then
                    write(*,*) "ERROR READING ATOM NUMBER",i
                endif
                write(*,*) "Read only position for atom",i
            endif
         endif
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
         write(6,*) 'possible choices are:'
         write(6,'(3XA)') type_simu
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
subroutine get_bool(io,fname,vname,val)
use m_vector
implicit none
logical, intent(inout) :: val
integer, intent(in) :: io
character(len=*), intent(in) :: vname,fname
!internal
logical :: success
integer :: stat
character(len=len(vname))   ::  tmp

Call set_pos_entry(io,fname,vname,success)
if(success)then
    read(io,*,iostat=stat) tmp, val
    if(stat/=0)then
        write(error_unit,'(/5A)') 'Failed to read "',vname,'" in file "',fname,'", but keyword is given'
        STOP "Fix input"
    endif
    write(output_unit,'(A,A30,3A,L)') 'Found entry for    ','"'//vname//'"',' in file "',fname,'", using value:         ', val
else
    write(output_unit,'(A,A30,3A,L)') 'No entry found for ','"'//vname//'"',' in file "',fname,'", using default value: ', val
endif
Call check_further_entry(io,fname,vname)
end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get a REAL number parameter (check the string and so on)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_real(io,fname,vname,val)
use m_vector
implicit none
real(kind=8), intent(inout) :: val
integer, intent(in) :: io
character(len=*), intent(in) :: vname,fname
!internal
logical :: success
integer :: stat
character(len=len(vname))   ::  tmp

Call set_pos_entry(io,fname,vname,success)
if(success)then
    read(io,*,iostat=stat) tmp, val
    if(stat/=0)then
        write(error_unit,'(/5A)') 'Failed to read "',vname,'" in file "',fname,'", but keyword is given'
        STOP "Fix input"
    endif
    write(output_unit,'(A,A30,3A,E16.8)') 'Found entry for    ','"'//vname//'"',' in file "',fname,'", using value:         ', val
else
    write(output_unit,'(A,A30,3A,E16.8)') 'No entry found for ','"'//vname//'"',' in file "',fname,'", using default value: ', val
endif
Call check_further_entry(io,fname,vname)
end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get a INTEGER number parameter (check the string and so on)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_int(io,fname,vname,val)
use m_vector
implicit none
integer, intent(inout)          :: val
integer, intent(in)             :: io
character(len=*), intent(in)    :: vname,fname
logical                     :: success
integer                     :: stat
character(len=len(vname))   ::  tmp

Call set_pos_entry(io,fname,vname,success)
if(success)then
    read(io,*,iostat=stat) tmp, val
    if(stat/=0)then
        write(error_unit,'(/5A)') 'Failed to read "',vname,'" in file "',fname,'", but keyword is given'
        STOP "Fix input"
    endif
    write(output_unit,'(A,A30,3A,I12)') 'Found entry for    ','"'//vname//'"',' in file "',fname,'", using value:         ', val
else
    write(output_unit,'(A,A30,3A,I12)') 'No entry found for ','"'//vname//'"',' in file "',fname,'", using default value: ', val
endif
Call check_further_entry(io,fname,vname)
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get a INTEGER number parameter (check the string and so on)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_int8(io,fname,vname,val)
use m_vector
implicit none
integer(8), intent(inout)       :: val 
integer, intent(in)             :: io
character(len=*), intent(in)    :: vname,fname
logical                     :: success
integer                     :: stat
character(len=len(vname))   ::  tmp

Call set_pos_entry(io,fname,vname,success)
if(success)then
    read(io,*,iostat=stat) tmp, val
    if(stat/=0)then
        write(error_unit,'(/5A)') 'Failed to read "',vname,'" in file "',fname,'", but keyword is given'
        STOP "Fix input"
    endif
    write(output_unit,'(A,A30,3A,I12)') 'Found entry for    ','"'//vname//'"',' in file "',fname,'", using value:         ', val
else
    write(output_unit,'(A,A30,3A,I12)') 'No entry found for ','"'//vname//'"',' in file "',fname,'", using default value: ', val
endif
Call check_further_entry(io,fname,vname)
end subroutine 



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

subroutine get_vec1d_real(io,fname,vname,val)
    real(8), intent(inout) :: val(:)
    integer, intent(in) :: io
    character(len=*), intent(in) :: vname,fname
    !internal
    logical :: success
    integer :: stat
    character(len=len(vname))   ::  tmp
   
    Call set_pos_entry(io,fname,vname,success)
    if(success)then
        read(io,*,iostat=stat) tmp, val
        if(stat/=0)then
            write(error_unit,'(/3A)') 'Failed to read ',vname,' but keyword is given'
            STOP "Fix input"
        endif
        if(size(val)<4)then
            write(output_unit,'(A,A30,3A,3E16.8)') 'Found entry for    ','"'//vname//'"',' in file "',fname,'", using value:         ', val
        else
            write(output_unit,'(A,A30,3A)') 'Found entry for    ','"'//vname//'"',' in file "',fname,'"'
            write(output_unit,'(3XA)')      '   The values are:'
            write(output_unit,'(18X,3E16.8)') val
        endif
    else
        if(size(val)<4)then
            write(output_unit,'(A,A30,3A,3E16.8)') 'No entry found for ','"'//vname//'"',' in file "',fname,'", using default value: ', val
        else
            write(output_unit,'(A,A30,3A)') 'No entry found for ','"'//vname//'"',' in file "',fname,'"'
            write(output_unit,'(3XA)')      '   Using the default values:'
            write(output_unit,'(18X,3E16.8)') val
        endif
    endif
    Call check_further_entry(io,fname,vname)
end subroutine 

subroutine get_vec1d_int(io,fname,vname,vec)
    integer, intent(inout) :: vec(:)
    integer, intent(in) :: io
    character(len=*), intent(in) :: vname,fname
    !internal
    logical :: success
    integer :: stat
    character(len=len(vname))   ::  tmp
   
    Call set_pos_entry(io,fname,vname,success)
    if(success)then
        read(io,*,iostat=stat) tmp, vec
        if(stat/=0)then
            write(error_unit,'(/3A)') 'Failed to read ',vname,' but keyword is given'
            STOP "Fix input"
        endif
    else
        write(output_unit,'(/2A)') 'No entry found for ',vname
        write(output_unit,'(A)') 'Using default value:'
        write(output_unit,'(3I10)') vec
    endif
    Call check_further_entry(io,fname,vname)
end subroutine 



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

function get_NB_columns(fname) result(ncol)
implicit none
character(len=*), intent(in) :: fname
integer :: ncol
! internal
integer :: npoint,io,i,nspace,fin
character(len=400) :: str
logical :: test
real(8),allocatable :: dumy(:)

ncol=0
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

  npoint=0
  do i=1,len(str)
    if (str(i:i) == '.') npoint=npoint+1
  enddo

  nspace=0
  do i=1,len(str)
    if (str(i:i) == ' ') nspace=nspace+1
  enddo

  ncol=npoint
  if (mod(nspace,npoint-1).eq.0) ncol=npoint

allocate(dumy(ncol),source=0.0d0)

rewind(io)
do
  read (io,*,iostat=fin) dumy
  if (fin .gt. 0) STOP 'ERROR in reading file in get_NB_columns'
  if (fin .lt. 0) exit
enddo
close(io)

write(output_unit,'(/a,I4,a/)') 'the file has ',ncol,' columns'

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
