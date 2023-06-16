module m_sym_dense
use m_symmetry_base
use m_sym_utils, only : look_translation
use m_derived_types, only : t_cell
use m_type_lattice, only : lattice
use m_lattice_position
use, intrinsic :: iso_fortran_env, only : output_unit
use m_invert, only : invert
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
class(pt_grp_dense), intent(in)    :: this
real(kind=8)       , intent(in)    :: areal(3,3)
integer            , intent(out)   :: number_sym
integer            , intent(inout) :: sym_index(:)
integer      ,intent(in)     :: dim_lat(:)
logical      ,intent(in)     :: periodic(:)

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

    found=look_translation(rtest,transpose(areal),(/.true.,.true.,.true./),dim_lat,tol_sym=this%tol_sym)

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
subroutine get_pt_sym(this,my_lattice,number_sym,sym_index,my_motif,sym_translation)
use m_vector, only : norm
class(pt_grp_dense), intent(in)   :: this
type(lattice), intent(in)         :: my_lattice
integer      , intent(inout)      :: number_sym
integer,intent(inout),allocatable :: sym_index(:)
type(t_cell) , intent(in)         :: my_motif
real(8),intent(inout),allocatable :: sym_translation(:,:)
!internal
integer :: natom,i,j,i_sim,k,nattype,size_pos,ii,n_at_under_type,j_atom
integer,allocatable ::  new_index(:),mask_index(:),ind(:)
real(8) :: test_vec(3),pos(3),sym_mat(3,3),test_pos(3)
integer,allocatable :: table_of_equivalence(:,:,:)
logical :: found

sym_mat=0.0d0
natom=size(my_motif%atomic)
allocate(new_index(number_sym),source=0)
allocate(mask_index(number_sym),source=0)
allocate(table_of_equivalence(3,natom,number_sym),source=0)

write(output_unit,'(I4,a)') natom, ' atoms found in the unit cell'

nattype=my_motif%n_attype    ! number of atom type

j_atom=0
do i=1,nattype  ! loop over the number of atom type

   pos=my_motif%atomic(i)%position
   call my_motif%ind_attype(i,ind)
   n_at_under_type=size(ind)

   if ((norm(pos).lt.this%tol_sym).and.(n_at_under_type.eq.1)) then
           mask_index=1       ! update the mask
       cycle    !the atom is at the origin so the symmetries are the symmetries of the lattice
   endif

   do ii=1,n_at_under_type  ! loop over all the atoms of the same nature

      pos=my_motif%atomic(ind(ii))%position ! take a test position
      j_atom=j_atom+1

      do j=1,number_sym    ! loop over all the symmetries
         i_sim=sym_index(j)
         if (i_sim.eq.0) cycle  ! if the symmetry does not apply to the lattice cycle
         table_of_equivalence(1,j_atom,i_sim)=i
         table_of_equivalence(2,j_atom,i_sim)=ii

         sym_mat=this%rotmat(i_sim)%mat
         test_vec=matmul(sym_mat,pos)    ! image of the test position by the symmetry operations

         do k=1,n_at_under_type    ! check if all atoms of the same type can be found within the symmetry operation
            test_pos=my_motif%atomic(ind(k))%position

            found=look_translation(test_vec,my_lattice%areal,my_lattice%periodic,my_lattice%dim_lat,test_pos,tol_sym=this%tol_sym)

            if (found) then
               mask_index(j)=mask_index(j)+1
               table_of_equivalence(3,j_atom,i_sim)=k
               exit
            endif
         enddo

      enddo
   enddo
   mask_index=mask_index/n_at_under_type
enddo

sym_index=sym_index*mask_index

j=0
do i=1,number_sym

   if (sym_index(i).eq.0) cycle
   j=j+1
   new_index(j)=sym_index(i)

enddo

number_sym=j
allocate(sym_translation(3,number_sym),source=0.0d0)
sym_index=new_index
write(6,'(a,I3,a)') 'space group has ',number_sym,'  symmetry operations'
write(6,*) (this%rotmat(sym_index(j))%name,j=1,number_sym)

! print the table of equivalence
do i=1,number_sym
   i_sim=sym_index(i)
   if (i_sim.eq.0) cycle
   write(6,'(/a)') '------------------------'
   write(6,'(2a)')'symmetry operation ',this%rotmat(i_sim)%name
   write(6,'(a)') '------------------------'
   do j=1,natom
      write(6,'(a,2I3,a,I3)')'atom type ',table_of_equivalence(1,j,i_sim), table_of_equivalence(2,j,i_sim), ' -> ', table_of_equivalence(3,j,i_sim)
   enddo
   write(6,'(a/)') '------------------------'
enddo

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

subroutine load(this,index_sym,all_sym,sym_translation)
 class(pt_grp_dense), intent(inout) :: this
 class(pt_grp)      , intent(in)    :: all_sym
 integer            , intent(in)    :: index_sym(:)
 real(8)            , intent(in)    :: sym_translation(:,:)

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
        write(io_sym,'(3(f14.10,3x))') this%rotmat(i)%translation
    enddo

    call close_file('symmetries.out',io_sym)

end subroutine

subroutine read_sym(this,fname)
    use m_io_files_utils
    class(pt_grp_dense),intent(inout) :: this
    character(len=*), intent(in) :: fname

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
        read(io_sym,*) this%rotmat(i)%translation
    enddo

    call close_file(fname,io_sym)

end subroutine

subroutine get_all_symetries(this,N,out)
    class(pt_grp_dense),intent(in)    :: this
    integer            ,intent(out)   :: N
    class(pt_grp)      ,intent(inout) :: out

    integer :: i,dim_table,j,k,fin
    type(pt_grp_dense) :: tmp_op
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
             if (all(abs(tmp_mat-tmp_op%rotmat(k)%mat).lt.this%tol_sym)) then
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
    class(pt_grp_dense),intent(in)    :: this
    class(pt_grp_dense),intent(inout) :: out

    type(pt_grp_dense) :: tmp_op
    integer :: dim_table,i

    dim_table=this%n_sym**2

    allocate(tmp_op%rotmat(dim_table))

    do i=1,this%n_sym
       tmp_op%rotmat(i)%mat=this%rotmat(i)%mat
       tmp_op%rotmat(i)%name=this%rotmat(i)%name
    enddo


end subroutine

end module
