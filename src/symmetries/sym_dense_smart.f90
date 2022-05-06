module m_sym_dense_smart
use m_symmetry_base
use m_sym_utils, only : look_translation
use m_derived_types, only : t_cell
use m_type_lattice, only : lattice
use m_lattice_position
use, intrinsic :: iso_fortran_env, only : output_unit
use m_invert, only : invert
implicit none

private
public :: pt_grp_dense_smart

type, extends(pt_grp) :: pt_grp_dense_smart
contains
   procedure :: init
   procedure :: load
   procedure :: apply_sym
   procedure :: get_N_sym
   procedure :: get_latt_sym
   procedure :: write_sym,read_sym
   procedure :: get_pt_sym
   procedure :: get_all_symetries

end type

contains

! find the symmetry of the unit cell
! Take the lattice vectors R(3,3) with basis vector ordered in lines
! apply a symmetry operation P on it P.R with vectors in column
! Calculate P.R-R
! find the linear combination so that (P.R-R) must be 0
! Continue until there are no symmetry operation left
!
subroutine get_latt_sym(this,areal,number_sym,sym_index,periodic,dim_lat)
class(pt_grp_dense_smart), intent(in) :: this
real(kind=8)       , intent(in)       :: areal(3,3)
integer            , intent(out)      :: number_sym
integer            , intent(inout)    :: sym_index(:)
integer            ,intent(in)        :: dim_lat(:)
logical            ,intent(in)        :: periodic(:)

integer :: i,j
real(kind=8) :: rtest(3,3)
logical :: found

number_sym=0
! loop over all symmetry
do i=1,size(sym_index)

!
! I multiply the rotation matrice written in line with the lattice vector written in lign too
! I have to take the transpose of alat to have it written in column
! the matrix rtest is then written in column in (column,line)
!
    rtest=matmul(this%rotmat(i)%mat,transpose(areal))

!
! I have to take the transpose of rtest to have it written in lign again
!

    found=look_translation(rtest,transpose(areal),periodic,dim_lat)

    if (found) then
       number_sym=number_sym+1
       sym_index(number_sym)=i
    endif

enddo

write(output_unit,'(/a,I4,a,I4)') 'number of lattice symmetries found  ', number_sym, '/',this%N_sym
do j=1,number_sym
   write(output_unit,'(a)') trim(this%rotmat(sym_index(j))%name)
   do i=1,3
      write(output_unit,'(3f8.3)') this%rotmat(sym_index(j))%mat(i,:)
   enddo
enddo

end subroutine get_latt_sym


! find the symmetry of each atomic position in the unit cell compatible with the lattice symmetry
! Take the lattice vectors R(3,3) with basis vector ordered in lines
! Take to atom position rho(3)
! Calculate R.rho
! apply a symmetry operation P on it P.(R.rho)
! Calculate P.(R.rho)-R
! find the linear combination so that (P.(R.rho)-R) must be 0
! Continue until there are no symmetry operation left
!
subroutine get_pt_sym(this,my_lattice,number_sym,sym_index,my_motif)
use m_vector, only : norm
class(pt_grp_dense_smart), intent(in)    :: this
type(lattice)            , intent(in)    :: my_lattice
integer                  , intent(inout) :: number_sym,sym_index(:)
type(t_cell)             , intent(in)    :: my_motif
!internal
integer :: natom,i,j,i_sim,k,nattype,size_pos,ii,npos,jj,max_sym,n_max_sym
integer,allocatable ::  new_index(:),mask_index(:,:),ind(:),ind_high_sym(:)
real(kind=8),allocatable :: high_sym_pos(:,:)
real(kind=8) :: test_vec(3),pos(3),sym_mat(3,3),test_pos(3),translation(3),origin(3)
logical :: found

sym_mat=0.0d0
natom=size(my_motif%atomic)
nattype=my_motif%n_attype
npos=1
do i=2,natom,1
   npos=npos*i
enddo

allocate(new_index(number_sym),source=0)
allocate(mask_index(number_sym,npos+6+1),source=0)
allocate(high_sym_pos(3,npos+6+1),source=0.0d0)

write(output_unit,'(I4,a)') natom, ' atoms found in the unit cell'


! find the high symmetry point
! 1. collect all atoms position
! 2. find the position of all possible high symmetry points (position of each atoms + inbetween each atoms) factorial(natom)
! 3. check all symmetry operations with possible translations

! 1.
high_sym_pos(:,1)=my_lattice%areal(1,:)/2.0d0
high_sym_pos(:,2)=my_lattice%areal(2,:)/2.0d0
high_sym_pos(:,3)=my_lattice%areal(3,:)/2.0d0
high_sym_pos(:,4)=(my_lattice%areal(1,:)+my_lattice%areal(2,:))/2.0d0
high_sym_pos(:,5)=(my_lattice%areal(2,:)+my_lattice%areal(3,:))/2.0d0
high_sym_pos(:,6)=(my_lattice%areal(1,:)+my_lattice%areal(3,:))/2.0d0

ii=6+1
do i=1,natom
   high_sym_pos(:,ii)=my_motif%atomic(i)%position
   ii=ii+1
enddo

!2.
do i=1,natom
   do j=i+1,natom
      high_sym_pos(:,ii)=(my_motif%atomic(i)%position+my_motif%atomic(j)%position)/2.0d0
      ii=ii+1
   enddo
enddo

!3.
do j=1,nattype  ! loop over the number of atom type
   call my_motif%ind_attype(j,ind)
   nattype=size(ind)

   do k=1,size(ind)  ! all atoms of the same type
      pos=my_motif%atomic(ind(k))%position ! take a test position

      do i=1,size(high_sym_pos,2)
         origin=high_sym_pos(:,i)   ! new origin of the unit cell
         translation=origin-pos     ! coordinate of the site on the new origin

         do ii=1,number_sym
            i_sim=sym_index(ii)
            if (i_sim.eq.0) cycle  ! if the symmetry does not apply to the lattice cycle
            sym_mat=this%rotmat(i_sim)%mat
            test_vec=matmul(sym_mat,translation)

            do jj=1,size(ind)
               test_pos=origin-my_motif%atomic(ind(jj))%position ! take a test position

               found=look_translation(test_vec,my_lattice%areal,(/.true.,.true.,.true./),my_lattice%dim_lat,test_pos)
               if (found) mask_index(ii,i)=mask_index(ii,i)+1
            enddo
         enddo

      enddo


   enddo
enddo

! treat the data
max_sym=0
do i=1,size(high_sym_pos,2)
   if(count(mask_index(:,i).eq.natom).ge.max_sym) max_sym=count((mask_index(:,i)).eq.natom)
enddo
write(output_unit,'(/a,I4)') 'maximum number of symmetry found: ',max_sym

j=0
do i=1,size(high_sym_pos,2)
   if(count(mask_index(:,i).eq.natom).eq.max_sym) j=j+1
enddo
n_max_sym=j
allocate(ind_high_sym(n_max_sym),source=0)
write(output_unit,'(a,I4,a)') 'found ',n_max_sym,' high symmetry positions'
write(output_unit,'(/a)') 'unit cell corresponding to the highest symmetries'
j=0
do i=1,size(high_sym_pos,2)
   if(count(mask_index(:,i).eq.natom).eq.max_sym) then
      j=j+1
      ind_high_sym(j)=i
      write(output_unit,'(/a,I4)') 'unit cell: ', j
      write(output_unit,'(a,3(E20.12E3,2x))') 'origin: ', high_sym_pos(:,i)
      write(output_unit,'(a)') 'atomic positions'
      do k=1,natom
         write(output_unit,'(a,3(E20.12E3,2x))') my_motif%atomic(k)%name,high_sym_pos(:,i)-my_motif%atomic(k)%position
      enddo
   endif
enddo

! check the unit cell - check of the origin should be changed
found=.false.
outer : do i=1,n_max_sym
   do k=1,natom
      if (norm(high_sym_pos(:,ind_high_sym(i))-my_motif%atomic(k)%position).lt.1.0d-8) then
         write(output_unit,'(a)') 'unit cell should not be changed'
         found=.true.
         exit outer
      endif
   enddo
enddo outer

if (.not.(found)) STOP 'restart the simulation with a better unit cell'

k=ind_high_sym(i)
mask_index=mask_index/natom        ! to get a series of 0 or 1
sym_index=sym_index*mask_index(:,k)

j=0
do i=1,number_sym

   if (sym_index(i).eq.0) cycle
   j=j+1
   new_index(j)=sym_index(i)

enddo

number_sym=j
sym_index=new_index
write(6,'(a,I3,a)') 'space group has ',number_sym,'  symmetry operations'
write(6,*) (this%rotmat(sym_index(j))%name,j=1,number_sym)

end subroutine get_pt_sym



subroutine init(this,N_sym)
 class(pt_grp_dense_smart), intent(out) :: this
 integer                  , intent(in)  :: N_sym

 integer :: i

 allocate(this%rotmat(N_sym))
 ! set to zero
 do i=1,N_sym
     this%rotmat(i)%mat = 0.d0
     this%rotmat(i)%name=""
 end do
 this%n_sym=N_sym

end subroutine

subroutine load(this,index_sym,all_sym)
 class(pt_grp_dense_smart), intent(inout) :: this
 class(pt_grp)            , intent(in)    :: all_sym
 integer                  , intent(in)    :: index_sym(:)

 integer :: i,j
 logical :: if_file

 if_file=.false.
 inquire(file='symmetries.in',exist=if_file)
! if file found
if (if_file) then
 write(output_unit,'(a/)') 'symmetries read from symmetries.out'
 call this%read_sym('symmetries.in')
else

! if file not found
 do i=1,this%N_sym
     j=index_sym(i)
     this%rotmat(i)%mat = transpose(all_sym%rotmat(j)%mat)
     this%rotmat(i)%name = all_sym%rotmat(j)%name
 end do
endif

end subroutine


function apply_sym(this,i,u) result(v)
 class(pt_grp_dense_smart), intent(in) :: this
 integer                  , intent(in) :: i
 real(8)                  , intent(in) :: u(3)
 real(8)                         :: v(3)

v=matmul(this%rotmat(i)%mat,u)

end function

function get_N_sym(this) result(N)
 class(pt_grp_dense_smart), intent(in) :: this
 integer                               :: N

 N=size(this%rotmat)

end function

!
! subroutine that write down the symmetry operations that were found
!
subroutine write_sym(this)
use m_io_files_utils
class(pt_grp_dense_smart), intent(in) :: this
! internal variables
integer :: io_sym
integer :: i,j

    io_sym=open_file_write('symmetries.out')

    write(io_sym,'(I4)') this%n_sym
    do i=1,this%n_sym
        write(io_sym,'(a)') '!'
        write(io_sym,'(a)') this%rotmat(i)%name
        do j=1,3
            write(io_sym,'(3(f14.10,3x))') this%rotmat(i)%mat(j,:)
        enddo
    enddo

    call close_file('symmetries.out',io_sym)

end subroutine

subroutine read_sym(this,fname)
    use m_io_files_utils
    class(pt_grp_dense_smart), intent(inout) :: this
    character(len=*)         , intent(in)    :: fname

    ! internal
    integer :: io_sym,i,j,N_sym


    io_sym=open_file_read(fname)

    read(io_sym,*)N_sym
    call this%init(N_sym)

    do i=1,n_sym
        read(io_sym,*)
        read(io_sym,*) this%rotmat(i)%name
        do j=1,3
            read(io_sym,*) this%rotmat(i)%mat(j,:)
        enddo
    enddo

    call close_file(fname,io_sym)

end subroutine

subroutine get_all_symetries(this,N,out)
    class(pt_grp_dense_smart),intent(in)    :: this
    integer                  ,intent(out)   :: N
    class(pt_grp)            ,intent(inout) :: out

    integer :: i,dim_table,j,k,fin
    type(pt_grp_dense_smart) :: tmp_op
    real(8) :: tmp_mat(3,3)
    logical :: found

    dim_table=this%n_sym**2
    allocate(tmp_op%rotmat(dim_table))

    do i=1,this%n_sym
       tmp_op%rotmat(i)%mat=this%rotmat(i)%mat
       tmp_op%rotmat(i)%name=this%rotmat(i)%name
    enddo

    fin=this%n_sym
    do i=1,this%n_sym
       do j=1,this%n_sym
       tmp_mat=matmul(tmp_op%rotmat(i)%mat,tmp_op%rotmat(j)%mat)

          found=.false.
          do k=1,fin
             if (all(abs(tmp_mat-tmp_op%rotmat(k)%mat).lt.1.0d-6)) then
                found=.true.
                exit
             endif
          enddo

          if (.not.found) then
             tmp_op%rotmat(fin+1)%mat=tmp_mat
             tmp_op%rotmat(fin+1)%name=trim(this%rotmat(i)%name)//'*'//trim(tmp_op%rotmat(j)%name)
             fin=fin+1
          endif
       enddo
    enddo

    write(output_unit,'(I5,2x,a)') fin,'different symmetry operations found'

    allocate(out%rotmat(fin))
    out%n_sym=fin
    N=fin
    do i=1,N
       out%rotmat(i)%mat=tmp_op%rotmat(i)%mat
       out%rotmat(i)%name=tmp_op%rotmat(i)%name
    enddo

end subroutine

!
!
! get the multiplication table for all the symmetry operations
!
subroutine get_multiplication_table(this,out)
    class(pt_grp_dense_smart),intent(in)    :: this
    class(pt_grp_dense_smart),intent(inout) :: out

    type(pt_grp_dense_smart) :: tmp_op
    integer :: dim_table,i

    dim_table=this%n_sym**2

    allocate(tmp_op%rotmat(dim_table))

    do i=1,this%n_sym
       tmp_op%rotmat(i)%mat=this%rotmat(i)%mat
       tmp_op%rotmat(i)%name=this%rotmat(i)%name
    enddo


end subroutine

end module
