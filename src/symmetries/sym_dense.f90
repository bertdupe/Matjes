module m_sym_dense
use m_symmetry_base
use m_sym_utils, only : look_translation
use m_derived_types, only : t_cell
use, intrinsic :: iso_fortran_env, only : output_unit
implicit none

private
public :: pt_grp_dense

type, extends(pt_grp) :: pt_grp_dense
contains
   procedure :: init
   procedure :: load
   procedure :: apply_sym
   procedure :: get_N_sym
   procedure :: get_latt_sym
   procedure :: write_sym,read_sym
   procedure :: get_pt_sym

end type

contains

! find the symmetry of the unit cell
! Take the lattice vectors R(3,3) with basis vector ordered in lines
! apply a symmetry operation P on it P.R with vectors in column
! Calculate P.R-R
! find the linear combination so that (P.R-R) must be 0
! Continue until there are no symmetry operation left
!
subroutine get_latt_sym(this,areal,number_sym,sym_index)
class(pt_grp_dense), intent(in)    :: this
real(kind=8)       , intent(in)    :: areal(3,3)
integer            , intent(inout) :: number_sym,sym_index(:)

integer :: i,j
real(kind=8) :: rtest(3,3)
logical :: found(3)

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

    do j=1,3
       found(j)=look_translation(rtest(:,j),areal,(/.true.,.true.,.true./),(/2,2,2/))
    enddo

    if (all(found)) then
       number_sym=number_sym+1
       sym_index(number_sym)=i
    endif

enddo

write(output_unit,'(/a,I4,a)') 'number of lattice symmetries found  ', number_sym, '/64'
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
subroutine get_pt_sym(this,areal,number_sym,sym_index,my_motif,periodic,dim_lat)
use m_vector, only : norm
class(pt_grp_dense), intent(in)     :: this
integer            , intent(inout)  :: number_sym,sym_index(:)
real(8)            , intent(in)     :: areal(3,3)
integer            , intent(in)     :: dim_lat(:)
logical            , intent(in)     :: periodic(:)
type(t_cell)       , intent(in)     :: my_motif
!internal
integer :: natom,i,j,i_sim,k
integer,allocatable ::  new_index(:),mask_index(:)
real(kind=8) :: test_vec(3),pos(3),sym_mat(3,3),test_pos(3)
logical :: found

sym_mat=0.0d0
natom=size(my_motif%atomic)
allocate(new_index(number_sym),source=0)
allocate(mask_index(number_sym),source=0)

write(output_unit,'(I4,a)') natom, ' atoms found in the unit cell'

! loop over all the symmetries
do j=1,number_sym
   i_sim=sym_index(j)
   if (i_sim.eq.0) cycle  ! if the symmetry does not apply to the lattice cycle

   sym_mat=this%rotmat(i_sim)%mat

   do i=1,natom  ! loop over the number of atoms

      pos=my_motif%atomic(i)%position
      test_vec=matmul(sym_mat,pos)

      do k=1,natom   ! check if any of the other atoms is equivalent

         test_pos=my_motif%atomic(k)%position
         found=look_translation(test_vec,areal,periodic,dim_lat,test_pos)

         if (found) then
             write(output_unit,'(3a,3f8.3)') "atom : ", trim(my_motif%atomic(i)%name), " at : ", my_motif%atomic(i)%position
             write(output_unit,'(3a,3f8.3)') "equivalent to atom : ", trim(my_motif%atomic(k)%name), " at : ", my_motif%atomic(k)%position
             write(output_unit,'(2a/)') "via symmetry : ", trim(this%rotmat(i_sim)%name)
             mask_index(j)=1
         endif

      enddo

   enddo

enddo

sym_index=sym_index*mask_index

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
 class(pt_grp_dense), intent(out) :: this
 integer            , intent(in)  :: N_sym

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
 class(pt_grp_dense), intent(inout) :: this
 class(pt_grp)      , intent(in)    :: all_sym
 integer            , intent(in)    :: index_sym(:)

 integer :: i,j

 do i=1,this%N_sym
     j=index_sym(i)
     this%rotmat(i)%mat = all_sym%rotmat(j)%mat
     this%rotmat(i)%name = all_sym%rotmat(j)%name
 end do

end subroutine


function apply_sym(this,i,u) result(v)
 class(pt_grp_dense), intent(in) :: this
 integer            , intent(in) :: i
 real(8)            , intent(in) :: u(3)
 real(8)                         :: v(3)

v=matmul(this%rotmat(i)%mat,u)

end function

function get_N_sym(this) result(N)
 class(pt_grp_dense), intent(in) :: this
 integer                         :: N

 N=size(this%rotmat)

end function

!
! subroutine that write down the symmetry operations that were found
!
subroutine write_sym(this)
use m_io_files_utils
class(pt_grp_dense), intent(in) :: this
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

subroutine read_sym(this)
    use m_io_files_utils
    class(pt_grp_dense),intent(inout) :: this

    ! internal
    integer :: io_sym,i,j,N_sym


    io_sym=open_file_read('symmetries.out')

    read(io_sym,*)N_sym
    call this%init(N_sym)

    do i=1,n_sym
        read(io_sym,*)
        read(io_sym,*) this%rotmat(i)%name
        do j=1,3
            read(io_sym,*) this%rotmat(i)%mat(j,:)
        enddo
    enddo

    call close_file('symmetries.out',io_sym)

end subroutine

end module
