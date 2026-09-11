module m_spglib
use m_symmetry_base
use m_lattice_position
use m_cell
use m_type_lattice
use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit, output_unit
#ifdef CPP_SPGLIB
use spglib_f08

implicit none

private
public :: sym_spglib

type, extends(pt_grp) :: sym_spglib
  type(SpglibDataset)        :: dset
  type(SpglibSpacegroupType) :: spgtype
  real(8)                    :: symprec
  integer                    :: max_num_sym=1000

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
class(sym_spglib)  , intent(in)       :: this
real(kind=8)       , intent(in)       :: areal(3,3)
integer            , intent(out)      :: number_sym
integer            , intent(inout)    :: sym_index(:)
integer            ,intent(in)        :: dim_lat(:)
logical            ,intent(in)        :: periodic(:)

number_sym = -1000

write(output_unit,'(a)') "routine not used in this setup"
continue

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
class(sym_spglib), intent(inout)         :: this
type(lattice)            , intent(in)    :: my_lattice
integer                  , intent(inout) :: number_sym
integer,intent(inout),allocatable        :: sym_index(:)
type(t_cell)             , intent(in)    :: my_motif
real(8),intent(inout),allocatable        :: sym_translation(:,:)
!internal

type(sym_spglib) :: test
real(8),allocatable :: positions(:)
integer, allocatable :: atom_types(:)
integer :: nat,nmag,i,int_number,intend
character(len=30) :: space

space = "                              "
nat=size(my_motif%atomic)
nmag=my_lattice%nmag
call my_motif%type_all_magnetic(atom_types)
call my_motif%pos_all_magnetic(positions)

test%dset = spg_get_dataset(lattice=my_lattice%areal, position=reshape(positions,(/3,nmag/)), types=atom_types, num_atom=nmag, symprec=this%tol_sym)

test%dset%spacegroup_number= spg_get_symmetry( rotation=test%dset%rotations, translation=test%dset%translations, max_size=test%max_num_sym, &
                            & lattice=my_lattice%areal,position=reshape(positions,(/3,nmag/)), types=atom_types, &
                            &  num_atom=nmag, symprec=this%tol_sym)

test%n_sym=size(test%dset%rotations,3)
allocate(test%rotmat(test%n_sym))
do i=1,test%n_sym
   test%rotmat(i)%mat=real(test%dset%rotations(:,:,i))     ! different order in C and fortran - transpose necessary
   test%rotmat(i)%translation=test%dset%translations(:,i)
   test%rotmat(i)%name=" "
enddo

int_number=spg_get_international( symbol=test%dset%international_symbol, lattice=my_lattice%areal, position=reshape(positions,(/3,nmag/)), &
                             & types=atom_types, num_atom=nmag, symprec=this%tol_sym)

if (int_number.eq.0) STOP "internal symbol could not be found in sym_spglib"

int_number=spg_get_schoenflies( symbol=test%spgtype%schoenflies, lattice=my_lattice%areal, &
                              & position=reshape(positions,(/3,nmag/)), types=atom_types, &
                              & num_atom=nmag, symprec=this%tol_sym );

if (int_number.eq.0) STOP "schoenflies group could not be found in sym_spglib"

intend=1
write(output_unit,'(a, "space_group: ", i3)') space(1:intend*2), test%dset%spacegroup_number
write(output_unit,'(a, "international: ", a, a)' ) space(1:intend*2), trim(test%dset%international_symbol)
write(output_unit,'(a, "schoenflies: ", a)') space(1:intend*2), trim(test%spgtype%schoenflies)

call write_syminfo( max_num_sym=test%max_num_sym, num_atom=nmag, lattice=my_lattice%areal, symprec=this%tol_sym, &
                              & atom_types=atom_types, positions=reshape(positions,(/3,nmag/)), &
                              & rotations=test%dset%rotations, translations=test%dset%translations, &
                              & international=test%dset%international_symbol, schoenflies=test%spgtype%schoenflies )

end subroutine get_pt_sym



subroutine init(this,N_sym)
 class(sym_spglib), intent(out) :: this
 integer          , intent(in)  :: N_sym

 integer :: i

 allocate(this%rotmat(N_sym))
 ! set to zero
 do i=1,N_sym
     this%rotmat(i)%mat = 0.0d0
     this%rotmat(i)%name=""
 end do
 this%n_sym=N_sym

end subroutine

subroutine load(this,index_sym,all_sym,sym_translation)
 class(sym_spglib), intent(inout) :: this
 class(pt_grp)    , intent(in)    :: all_sym
 integer          , intent(in)    :: index_sym(:)
 real(8)          , intent(in)    :: sym_translation(:,:)

continue

end subroutine


function apply_sym(this,i,u) result(v)
 class(sym_spglib), intent(in) :: this
 integer          , intent(in) :: i
 real(8)          , intent(in) :: u(3)
 real(8)                       :: v(3)

v=matmul(this%rotmat(i)%mat,u)+this%rotmat(i)%translation

end function

function get_N_sym(this) result(N)
 class(sym_spglib), intent(in) :: this
 integer                       :: N

 N=size(this%rotmat)

end function

!
! subroutine that write down the symmetry operations that were found
!
subroutine write_sym(this)
use m_io_files_utils
class(sym_spglib), intent(in) :: this
! internal variables
integer :: io_sym
integer :: i,j

    io_sym=open_file_write('symmetries.out')

    write(io_sym,'(I10)') this%n_sym
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
    class(sym_spglib), intent(inout) :: this
    character(len=*) , intent(in)    :: fname

    ! internal
    integer :: io_sym,i,j,N_sym


    io_sym=open_file_read(fname)

    read(io_sym,*) N_sym
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
    class(sym_spglib),intent(in)    :: this
    integer          ,intent(out)   :: N
    class(pt_grp)    ,intent(inout) :: out

    N=-1000

    continue

end subroutine

!
!
! get the multiplication table for all the symmetry operations
!
subroutine get_multiplication_table(this,out)
    class(sym_spglib),intent(in)    :: this
    class(sym_spglib),intent(inout) :: out

    type(sym_spglib) :: tmp_op
    integer :: dim_table,i

    dim_table=this%n_sym**2

    allocate(tmp_op%rotmat(dim_table))

    do i=1,this%n_sym
       tmp_op%rotmat(i)%mat=this%rotmat(i)%mat
       tmp_op%rotmat(i)%name=this%rotmat(i)%name
    enddo


end subroutine


subroutine write_syminfo( max_num_sym, num_atom, &
     lattice, symprec, atom_types, positions, rotations, translations, &
     international, schoenflies )

  ! Arguments ------------------------------------
  ! scalars
  integer, intent(in) :: num_atom, max_num_sym
  real(8), intent(in) :: symprec
  ! arrays
  integer, intent(in), dimension(:) :: atom_types
  real(8), intent(in), dimension(:, :) :: lattice
  real(8), intent(in), dimension(:, :) :: positions
  integer, intent(in), dimension(:, :, :) :: rotations
  real(8), intent(in), dimension(:, :) :: translations
  character(len=*) :: international, schoenflies
  ! Local variables-------------------------------
  ! scalars
  integer :: i, j, counter, indent
  character(len=30) :: space
  ! arrays
  real(8), dimension(3, 3) :: lattice_t
  !**************************************************************************

  space = "                              "


  ! transpose due to array order difference between C and fortran
  lattice_t = transpose( lattice )

  write (output_unit,'(a, "atom-type:")') space(1:indent*2)
  do i = 1, num_atom
     write(output_unit,'(a, "- { type: ", i3, " }")') space(1:indent*2), atom_types(i)
  end do
  write(output_unit,'(a, "real-basis:")') space(1:indent*2)
  do i = 1, 3
     write(output_unit,'(a, "- [", f19.14, ", ", f19.14, ", ", f19.14, " ]")') space(1:indent*2), lattice(:, i)
  end do
  print('(a, "position:")'), space(1:indent*2)
  do i = 1, num_atom
     write(output_unit,'(a, "- [", f17.14, ", ", f17.14, ", ", f17.14, " ]")') space(1:indent*2), positions(:, i)
  end do
  write(output_unit,'(a, "operation:")') space(1:indent*2)
  do i = 1, size(rotations,3)
     write(output_unit,'(a, "- rotation: #", i4)') space(1:indent*2), i
     do j = 1, 3
        write(output_unit,'(a, "  - [", i3,",", i3,",", i3," ]")') space(1:indent*2), rotations(j,:, i)
     end do
     write(output_unit,'(a, "  translation: [ ", f10.7,", ", f10.7,", ", f10.7," ]")') space(1:indent*2), translations(:,i)
  end do

!
! part that does the reciprocal mesh
!

!  print('(a, "reciprocal-mesh: [", i3, ", ", i3, ", ", i3, " ]")'), space(1:indent*2), mesh(:)
!  print('(a, "- is_shift: [", i3, ", ", i3, ", ", i3, " ]")'), space(1:indent*2), is_shift(:)
!  print('(a, "- is_time_reversal: ", i3)'), space(1:indent*2), is_time_reversal
!
!  call spg_get_ir_reciprocal_mesh( num_ir_grid, grid_point, map, &
!       mesh, is_shift, is_time_reversal, lattice, positions, &
!       atom_types, num_atom, symprec )
!
!  print('(a, "- num-ir-grid-point:", i4)'), space(1:indent*2), num_ir_grid
!  print('(a, "- ir-grid-point:")'), space(1:indent*2)
!  counter = 0
!  sum_weight = 0
!  do i = 1, mesh(1)*mesh(2)*mesh(3)
!     if ( i == map(i) ) then
!        ! Ad-hoc and intuitive implementation of weight
!        weight = 0
!        do j = 1, mesh(1)*mesh(2)*mesh(3)
!           if ( i == map(j) ) then
!              weight = weight + 1
!           end if
!        end do
!
!        counter = counter + 1
!        sum_weight = sum_weight + weight
!
!        print('(a, "  - #", i4)'), space(1:indent*2), counter
!        print('(a, "    address: [", i3, ", ", i3, ", ", i3, " ]")'), space(1:indent*2), grid_point(:, map(i))
!        print('(a, "    weight: ", i4)'), space(1:indent*2), weight
!     end if
!  end do
  ! print('(a, "  sum_weight: ", i4)'), space(1:indent*2), sum_weight

end subroutine write_syminfo


#endif

end module m_spglib

