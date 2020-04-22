module m_setup_DM
use m_grp_sym
use m_derived_types, only : cell,lattice
use m_basic_types, only : symop
use m_vector, only : cross,norm


private
public :: setup_DM_vector, get_number_DMI
contains

!
! This is a bunch of subroutines that are absolutely horrible. The point is to calculate the DMI
! Do not look into this unless you can not sleep at night or you wish to have a migraine
!

integer function get_number_DMI(fname)
use m_io_files_utils
use m_io_utils
implicit none
character(len=*) :: fname
! internal variables
integer :: io
real(kind=8) :: val

io=open_file_read(fname)
get_number_DMI=count_variables(io,'DMI_',fname)

if (get_number_DMI.gt.0) then
  call get_parameter(io,fname,'DMI_1',val)
  if (abs(val).lt.1.0d-8) get_number_DMI=0
endif
call close_file(fname,io)

end function




!
! Calculate the possible DM candidates
!

subroutine setup_DM_vector(indexNN,n_DMI,my_lattice,my_motif,DM_vector,tabledist)
implicit none
integer, intent(in) :: indexNN(:,:),n_DMI
type(cell), intent(in) :: my_motif
type(lattice), intent(in) :: my_lattice
real(kind=8), intent(in) :: tabledist(:,:)
real(kind=8), intent(inout) :: DM_vector(:,:,:)
! internal
integer :: phase,n_sym,atom_all_shells,k,j,i
logical :: inquire_file
type(symop), allocatable :: symmetries(:)

phase=size(DM_vector,3)

inquire_file=.false.
inquire(file='symmetries.out',exist=inquire_file)
if (.not.inquire_file) then
   call get_group(my_lattice%areal,my_motif,my_lattice%boundary,my_lattice%dim_lat)
endif

n_sym=get_num_sym_file()
allocate(symmetries(n_sym))
call read_symmetries(n_sym,symmetries)

j=1
do i=1,n_DMI
  atom_all_shells=indexNN(i,1)
  DM_vector(j:j+atom_all_shells-1,:,:)=get_DM_vectors(i,atom_all_shells,my_lattice%areal,my_motif,n_sym,symmetries,tabledist(i,1),my_lattice%boundary)
  j=j+atom_all_shells
enddo

write(6,'(a)') ''
write(6,'(a,2x,I4,2x,a)') 'find', atom_all_shells ,'DM vectors'
do j=1,phase
   do k=1,size(DM_vector,1)
      write(6,'(3(f8.5,2x))') DM_vector(k,:,j)
   enddo
enddo
write(6,'(a)') ''

end subroutine setup_DM_vector



!
! Calculate all the DMI in the unit cell
!
function get_DM_vectors(shell_number,atom_all_shells,areal,my_motif,n_sym,symmetries,distance,periodic)
implicit none
integer, intent(in) :: atom_all_shells,n_sym,shell_number
real(kind=8), intent(in) :: areal(:,:),distance
type(cell), intent(in) :: my_motif
type(symop), intent(in) :: symmetries(:)
real(kind=8) :: get_DM_vectors(atom_all_shells,3,1)
logical, intent(in) :: periodic(:)
! internal
integer :: i,natom,n_at_mag,n_at_nonmag,j,k,l,n
integer :: n_at_nonmag_incell,n_at_mag_inshell
real(kind=8) :: R1(3),R2(3),DMV(3),test,test_pos_A(3),test_pos_B(3),dumy
real(kind=8),allocatable :: R_mag(:,:,:),R_non_mag(:,:,:)
real(kind=8),allocatable :: pos_mag_atom_in_shell(:,:),pos_nonmag_atom_in_shell(:,:)
logical :: if_present,exists

natom=size(my_motif%atomic)
exists=.false.

n_at_mag=count(my_motif%i_mom(:),1)
n_at_nonmag=natom-n_at_mag
get_DM_vectors=0.0d0

!
! find atom in the shell 0
!
allocate(R_mag(3,n_at_mag,2),R_non_mag(3,n_at_nonmag,2))
R_mag=0.0d0
R_non_mag=0.0d0
call find_mag_and_nonmag_atoms(0,my_motif,areal,0.0d0,periodic,R_mag(:,:,1),R_non_mag(:,:,1))

call find_mag_and_nonmag_atoms(shell_number,my_motif,areal,distance,periodic,R_mag(:,:,2),R_non_mag(:,:,2))

! find all magnetic atoms in the unit cell
! the number of directions for the different atoms are put in pos_mag_atom_in_shell
! for that apply symmetry operations
n_at_mag_inshell=0
do j=1,n_at_mag
  do i=1,n_at_mag
    n_at_mag_inshell=n_at_mag_inshell+find_all_equivalent_atom(R_mag(:,i,2)-R_mag(:,j,1),n_sym,symmetries)
  enddo
enddo
allocate(pos_mag_atom_in_shell(3,n_at_mag_inshell))
call find_all_atom(pos_mag_atom_in_shell,R_mag(:,:,2),R_mag(:,:,1),n_sym,symmetries,areal)

! find all equivalent non magnetic atoms
! for that apply symmetry operations
n_at_nonmag_incell=0
do j=1,n_at_mag
  do i=1,n_at_nonmag
    n_at_nonmag_incell=n_at_nonmag_incell+find_all_equivalent_atom(R_non_mag(:,i,2)-R_mag(:,j,2),n_sym,symmetries)
  enddo
enddo

allocate(pos_nonmag_atom_in_shell(3,n_at_nonmag_incell))
call find_all_atom(pos_nonmag_atom_in_shell,R_non_mag(:,:,2),R_mag(:,:,2),n_sym,symmetries,areal)

! for all magnetic atoms in the unit cell
! find the unique non-magnetic atom that verifies the Moriya rules
l=1
do i=1,n_at_mag_inshell

   DMV=check_moryia_rules(pos_mag_atom_in_shell(:,i),pos_nonmag_atom_in_shell)

   if_present=.false.
   do j=1,l
      if (norm(get_DM_vectors(l,:,1)-DMV).lt.1.0d-8) if_present=.true.
   enddo

   if (.not.(if_present)) then
      get_DM_vectors(l,:,1)=DMV
      l=l+1
   endif
enddo

inquire(file='DM-2donly',exist=exists)
if (exists) then
   get_DM_vectors(:,3,:)=0.0d0
   do j=1,size(get_DM_vectors,3)
      do i=1,size(get_DM_vectors,1)
         dumy=norm(get_DM_vectors(i,:,j))
         get_DM_vectors(i,:,j)=get_DM_vectors(i,:,j)/dumy
      enddo
   enddo
endif

end function get_DM_vectors

















subroutine find_mag_and_nonmag_atoms(shell_number,my_motif,areal,distance,periodic,R_mag,R_non_mag)
implicit none
integer, intent(in) :: shell_number
real(kind=8), intent(in) :: areal(:,:),distance
type(cell), intent(in) :: my_motif
logical, intent(in) :: periodic(:)
real(kind=8), intent(inout) :: R_mag(:,:),R_non_mag(:,:)
! internal
integer :: i,j,k,natom,i1,i2,i3
real(kind=8) :: test(3)

natom=size(my_motif%atomic)

!
! find a set a coordinates to bring the atoms in the shell number shell_number
!
i=0
j=0
k=0
do i1=0,shell_number
  do i2=0,shell_number
    do i3=0,shell_number
       if (periodic(1)) i=i1
       if (periodic(2)) j=i2
       if (periodic(3)) k=i3
       test=matmul(transpose(areal),real((/i,j,k/)))
       if (abs( norm(test) -distance).lt.1.0d-8) exit
     enddo
     if (abs( norm(test) -distance).lt.1.0d-8) exit
  enddo
  if (abs( norm(test) -distance).lt.1.0d-8) exit
enddo

j=0
k=0
do i=1,natom
   if (abs(my_motif%atomic(i)%moment).gt.1.0d-8) then
      j=j+1
      R_mag(:,j)=matmul(transpose(areal),real((/i1,i2,i3/)))+my_motif%atomic(i)%position(:)
   else
      k=k+1
      R_non_mag(:,k)=matmul(transpose(areal),real((/i1,i2,i3/)))+my_motif%atomic(i)%position(:)
   endif

enddo

end subroutine find_mag_and_nonmag_atoms
!
! function that finds ne number of equivalent atoms in the unit cell
!

function find_all_equivalent_atom(R,n_sym,symmetries)
implicit none
real(kind=8), intent(in) :: R(:)
integer, intent(in) :: n_sym
type(symop), intent(in) :: symmetries(:)
integer :: find_all_equivalent_atom
! internal variables
integer :: i,j,l
real(kind=8) :: test_pos(3)
logical :: if_present
real(kind=8),allocatable :: atom_position(:,:)

allocate(atom_position(3,n_sym+1))
atom_position=0.0d0
test_pos=0.0d0

j=2
atom_position(:,1)=R
do i=1,n_sym
   test_pos=matmul(symmetries(i)%mat,R(:))

   if_present=.false.

   do l=1,j-1
      if (norm(atom_position(:,l)-test_pos).lt.1.0d-8) if_present=.true.
   enddo

   if (.not.(if_present)) then
      atom_position(:,j)=test_pos
      j=j+1
   endif

enddo

find_all_equivalent_atom=j-1

end function find_all_equivalent_atom

!
! function that finds all the equivalent atom in the unitcell
! choose an atom type (magnetic or not)
! calculate the position of all the atoms which at a distance R_1 of one atom.
!

subroutine find_all_atom(pos_atom_in_shell,R,Ref,n_sym,symmetries,areal)
implicit none
real(kind=8), intent(inout) :: pos_atom_in_shell(:,:)
real(kind=8), intent(in) :: R(:,:),Ref(:,:),areal(:,:)
integer, intent(in) :: n_sym
type(symop), intent(in) :: symmetries(:)
! internal
real(kind=8) :: test_pos(3)
integer :: j,l,k,i,n
logical :: if_present

j=2
pos_atom_in_shell(:,1)=R(:,1)-Ref(:,1)

do n=1,size(Ref,2)
  do i=1,size(R,2)

    do k=1,n_sym

      test_pos=matmul(symmetries(k)%mat,R(:,i)-Ref(:,n))

      if_present=.false.
      do l=1,j-1
        if (norm(pos_atom_in_shell(:,l)-test_pos).lt.1.0d-8) if_present=.true.
      enddo

      if (.not.(if_present)) then
        pos_atom_in_shell(:,j)=test_pos
        j=j+1
      endif

      if ((j-1).gt.size(pos_atom_in_shell,2)) stop 'error in setup_DM - find_all_atom'
    enddo
  enddo
enddo

write(6,'(a)') ''
write(6,'(a)') 'for atom at position:'
write(6,'(3(2x,f8.5))') pos_atom_in_shell(:,1)
write(6,'(a,2x,I4,2x,a)') 'find', size(pos_atom_in_shell,2)-1, 'equivalent atom in the unit cell'
do i=2,size(pos_atom_in_shell,2)
  write(6,'(3(2x,f8.5))') pos_atom_in_shell(:,i)
enddo
write(6,'(a)') ''

end subroutine find_all_atom

!
! function that checks if the DM vector verifies the Moriya rules
! choose a non magnetic direction (AB) and a non-magnetic atom direction at position C
! when the projection of the non-magnetic direction on AB (called C) is:
! (1) - a center of inversion for the triangle ABC then D=0
! (2) - when a mirror plane perpendicular to AB passes though C that DMI is in this plane
! (3) - when there is a mirror plane including A and B then DMI perpendicular to the mirror plane
! (4) - when there is a two-fold rotation axis perpendicular to AB and going through C and D perpendicular to the two fold axis
! (5) - when there is a n-fold rotation axis (n>1) along AB then D parallel to AB
!
function check_moryia_rules(R,R_nonmag)
implicit none
real(kind=8), intent(in) :: R(:),R_nonmag(:,:)
real(kind=8) :: check_moryia_rules(3)
! internal
real(kind=8) :: DM_candidate(3),test
integer :: i

check_moryia_rules=0.0d0
DM_candidate=0.0d0

do i=1,size(R_nonmag,2)

! rule (1) collinearity and inversion of symmetry
   if ((norm(cross(R,R_nonmag(:,i),1,3)).lt.1.0d-8).and.(norm(R_nonmag(:,i)-R/2.0d0).lt.1.0d-8)) cycle

! rule (2) and (3) DMI must be perpendicular to R. By construction it is always true
   DM_candidate=cross(R,R_nonmag(:,i),1,3)

   check_moryia_rules=check_moryia_rules+DM_candidate

enddo

test=norm(check_moryia_rules)
if (test.lt.1.0d-8) stop 'cannot find DM vector - check_moryia_rules '

check_moryia_rules=check_moryia_rules/test

end function

end module
