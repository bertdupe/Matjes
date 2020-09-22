module m_init_modes
use m_derived_types
use m_io_files_utils

interface get_init_modes
  module procedure init_3Dmodes
end interface get_init_modes

interface get_fast_index
!  module procedure get_fast_index_1D
!  module procedure get_fast_index_2D
  module procedure get_fast_index_3D
end interface get_fast_index

private
public :: get_init_modes

contains

!------------------------------------
! find the fast and the slow variable to make sure that the magetis texture is read correctly
!------------------------------------
subroutine get_fast_index_3D(fname,u,Ilat,N,lattice_vec,dim_lat,n_column)
implicit none
type(int_pointer),intent(out) :: u(:),N(:)
integer, target, intent(inout) :: Ilat(:),dim_lat(:)
integer, intent(in) :: n_column
real(kind=8), intent(in) :: lattice_vec(:,:)
character(len=*), intent(in) :: fname
! internal variables
real(kind=8) :: position(n_column,dim_lat(1),dim_lat(2),dim_lat(3))
integer :: i_x,i_y,i_z,j_lat,pos_origin(3)
integer :: i,i_1,i_2,i_3,order,npos,nsize
integer :: io
real(kind=8) :: dist_dum1,dist_dum2,origin(n_column),vec_test(n_column)
logical :: check

check=.False.
nsize=size(u)
Ilat(:)=1

do i=1,nsize
  u(i)%p=>Ilat(i)
  N(i)%p=>dim_lat(i)
enddo

dist_dum1=10.0d0
dist_dum2=0.0d0
origin=0.0d0
vec_test=0.0d0
order=0

io=open_file_read(fname)

!
! if fname does not exists, return with the defaut configuration in line 42 to 45
!
if (io.lt.0) then
  write(6,'(10a)') fname,' does not exist'
  return
endif

! read the position in cartesian coordinates
! the default is x,y,z in this order
do i_z=1,dim_lat(3)
   do i_y=1,dim_lat(2)
      do i_x=1,dim_lat(1)

        read(io,*) (position(j_lat,i_x,i_y,i_z),j_lat=1,n_column)

! find the origin of the coordinate system
        dist_dum2=sqrt(sum(position(:,i_x,i_y,i_z))**2)
        if (dist_dum2.lt.dist_dum1) then
           dist_dum1=dist_dum2
           origin=position(:,i_x,i_y,i_z)
           pos_origin(1)=i_x
           pos_origin(2)=i_y
           pos_origin(3)=i_z
        endif

      enddo
   enddo
enddo

! try the vector in the first direction

i_x=pos_origin(1)
i_y=pos_origin(2)
i_z=pos_origin(3)

do i_3=0,1
   if (i_3+i_z.gt.dim_lat(3)) cycle
   do i_2=0,1
      if (i_2+i_y.gt.dim_lat(2)) cycle
      do i_1=0,1
      if (i_1+i_x.gt.dim_lat(1)) cycle

      vec_test=position(:,i_x+i_1,i_y+i_2,i_z+i_3)-position(:,i_x,i_y,i_z)
       do i=1,3
         dist_dum1=sqrt(sum((vec_test-lattice_vec(i,:))**2))

         if (dist_dum1.lt.1.0d-8) then
           write(6,'(a)') 'lattice vector found'
           order=order+1
           write(6,'(a,I4,a)') 'vector number ',order,' found'
           npos=i_1*1+i_2*2+i_3*3

           write(6,'(2(a,I3,x),2a)') 'the index  ',i,' is the coordinate  ',npos,' in file  ', fname
           u(i)%p=>Ilat(npos)
           N(i)%p=>dim_lat(npos)

         endif

       enddo

      enddo
   enddo
enddo

rewind(io)
do i_3=1,N(3)%p
   u(3)%p=i_3
   do i_2=1,N(2)%p
      u(2)%p=i_2
      do i_1=1,N(1)%p
        u(1)%p=i_1

        i_x=Ilat(1)
        i_y=Ilat(2)
        i_z=Ilat(3)

      read(io,*) (position(j_lat,i_x,i_y,i_z),j_lat=1,n_column)

      enddo
   enddo
enddo

call close_file(fname,io)

end subroutine get_fast_index_3D

!------------------------------------
! find the column corresponding to the positions and the spins
!------------------------------------

subroutine get_col_spins(fname,Nstart,Nend,n_column,n_line)
use m_io_utils
implicit none
character(len=*), intent(in) :: fname
integer, intent(in) :: n_column,n_line
integer,intent(out) :: Nstart,Nend
!  internal
integer :: io,i,j_lat,Nstart_loc,Nend_loc
real(kind=8) :: dummy,diff
real(kind=8) :: test(n_column,n_line)

write(6,'(a)') 'locate the columns that contains the spins'

diff=1.0d-5
Nstart_loc=n_column-2
Nend_loc=n_column

io=open_file_read(fname)

do i=1,n_line
   read(io,*) (test(j_lat,i),j_lat=1,n_column)

   dummy=sum(test(Nstart_loc:Nend_loc,i)**2)

   if (abs(1.0d0-dummy).gt.diff) then
     Nstart_loc=Nstart_loc-1
     Nend_loc=Nend_loc-1
   endif
enddo

call close_file(fname,io)

write(6,'(a,I3,a,I3)') 'The  magnetic moments are located between column ',Nstart_loc,' and ',Nend_loc

Nstart=Nstart_loc
Nend=Nend_loc

end subroutine get_col_spins


subroutine init_3Dmodes(fname,my_lattice,motif)
use m_io_utils
use m_get_position
implicit none
type (lattice), intent(inout) :: my_lattice
type (cell), intent(in) :: motif
character(len=*), intent(in) :: fname
! internal variables
integer :: dim_lat(3),n_column,io,nmag,n_column_pos,io_pos
! pointer variable for internal use
type(int_pointer) :: u(3),N(3)
! slope variables
integer,target :: i_x,i_y,i_z,i_m,Ilat(4),Nstart,Nend
integer :: j_lat,i,i_1,i_2,i_3
real(kind=8) :: lattice_vec(3,3)
real(kind=8),allocatable :: dumy(:)
! local lattice
real(kind=8), allocatable :: local_lattice(:,:,:,:,:)

dim_lat=my_lattice%dim_lat
n_column=get_cols(fname)
nmag=count(motif%atomic(:)%moment.gt.0.0d0)
lattice_vec=my_lattice%areal
Ilat=0

do i=1,3
 nullify(u(i)%p,N(i)%p)
enddo

allocate(local_lattice(my_lattice%dim_mode,dim_lat(1),dim_lat(2),dim_lat(3),nmag))

if (n_column.eq.my_lattice%dim_mode) then

  ! read only the file containing the positions and assume that the magnetic moments are written
  ! with the same order
  n_column_pos=get_cols('positions.dat')
  call get_fast_index('positions.dat',u,Ilat,N,lattice_vec,dim_lat,n_column_pos)

elseif (n_column.gt.my_lattice%dim_mode+1) then

  ! read the init.config file and check the order of reading
  call get_col_spins(fname,Nstart,Nend,n_column,product(dim_lat)*nmag)
  ! read the init.config file and check the order of reading
  call get_fast_index(fname,u,Ilat,N,lattice_vec,dim_lat,Nstart-1)
else

  write(6,'(a)') 'WARNING - user defined format - I hope wou what you are doing'
  write(6,'(a)') 'check IO-modes line 268 for more information'

endif

n_column=get_cols(fname)
! Read the configurations
io=open_file_read(fname)

if (n_column.eq.my_lattice%dim_mode) then


do i_m=1,nmag
   do i_3=1,N(3)%p
      u(3)%p=i_3
      do i_2=1,N(2)%p
         u(2)%p=i_2
         do i_1=1,N(1)%p
           u(1)%p=i_1

           i_x=Ilat(1)
           i_y=Ilat(2)
           i_z=Ilat(3)
           read(io,*) (local_lattice(j_lat,i_x,i_y,i_z,i_m),j_lat=1,n_column)
         enddo
      enddo
   enddo
enddo

elseif (n_column.gt.my_lattice%dim_mode+1) then

allocate(dumy(Nstart-1))

do i_m=1,nmag
   do i_z=1,dim_lat(3)
      do i_y=1,dim_lat(2)
         do i_x=1,dim_lat(1)
           read(io,*) (dumy(i),i=1,Nstart-1),(local_lattice(j_lat,i_x,i_y,i_z,i_m),j_lat=1,3)
         enddo
      enddo
   enddo
enddo

deallocate(dumy)

else
!
! format to change by hand
!
!io_pos=open_file_read('positions.dat')
allocate(dumy(7))
local_lattice(4:5,:,:,:,:)=0.0d0
local_lattice(6,:,:,:,:)=20.0d0
local_lattice(1:2,:,:,:,:)=0.0d0
local_lattice(3,:,:,:,:)=1.0d0
do i_m=1,nmag
   do i_x=1,100
      do i_y=1,100
         do i_z=1,1
!           if ((i_x.le.100).and.(i_y.le.100)) then
             read(io,*) (dumy(i),i=1,3),(local_lattice(j_lat,i_x,i_y,i_z,i_m),j_lat=1,3),dumy(1)
!           else
!             read(io,*) (dumy(i),i=1,7)
!           endif

         enddo
      enddo
   enddo
enddo
deallocate(dumy)
endif

call close_file(fname,io)

! put the local_lattice variables in the lattice modes

do i_m=1,nmag
   do i_z=1,dim_lat(3)
      do i_y=1,dim_lat(2)
         do i_x=1,dim_lat(1)
        my_lattice%ordpar%l_modes(i_x,i_y,i_z,i_m)%w(:)=local_lattice(:,i_x,i_y,i_z,i_m)
         enddo
      enddo
   enddo
enddo

deallocate(local_lattice)

end subroutine init_3Dmodes

end module m_init_modes
