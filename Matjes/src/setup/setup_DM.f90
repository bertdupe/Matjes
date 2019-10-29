module m_setup_DM
use m_grp_sym
use m_derived_types, only : cell,lattice
use m_basic_types, only : symop
use m_vector, only : cross,norm

interface setup_DM
    module procedure setup_DM_1D
    module procedure setup_DM_2D
    module procedure setup_DM_3D
end interface setup_DM

private
public :: setup_DM_vector
contains

subroutine setup_DM_vector(indexNN,i,my_lattice,my_motif,DM_vector)
implicit none
integer, intent(in) :: indexNN(:,:),i
type (cell), intent(in) :: my_motif
type (lattice), intent(in) :: my_lattice
real(kind=8), intent(inout) :: DM_vector(:,:,:)
! internal
integer :: phase,n_sym,atom_all_shells,k,j
logical :: inquire_file
type(symop), allocatable :: symmetries(:)

phase=size(DM_vector,3)

inquire_file=.false.
inquire(file='symmetries.out',exist=inquire_file)
if (.not.inquire_file) then
   call get_group(my_lattice%areal,my_motif,my_lattice%boundary)
endif

n_sym=get_num_sym_file()
allocate(symmetries(n_sym))
call read_symmetries(n_sym,symmetries)

atom_all_shells=sum(indexNN(1:i,1))
if (.not.my_lattice%boundary(3)) DM_vector=get_DM_vectors(atom_all_shells,my_lattice%areal,my_motif,n_sym,symmetries)

#ifdef CPP_DEBUG
do j=1,phase
   do k=1,atom_all_shells
      write(6,'(3f7.5)') DM_vector(k,:,j)
   enddo
enddo
#endif

!if (size(my_lattice%world).eq.1) then
!    DM_vector=setup_DM(sum(indexNN(1:i,1)),my_lattice%areal,my_motif,phase)
!
!elseif (size(my_lattice%world).eq.2) then
!    DM_vector=setup_DM(sum(indexNN(1:i,1)),i,indexNN(:,1),my_lattice%areal,my_motif,my_lattice%world,phase)
!
!else
!    DM_vector=setup_DM(sum(indexNN(1:i,1)),my_lattice%areal,my_motif,phase)
!
!endif

end subroutine setup_DM_vector

!
! Calculate all the DMI in the uni cell
!
function get_DM_vectors(atom_all_shells,areal,my_motif,n_sym,symmetries)
implicit none
integer, intent(in) :: atom_all_shells,n_sym
real(kind=8), intent(in) :: areal(:,:)
type(cell), intent(in) :: my_motif
type(symop), intent(in) :: symmetries(:)
real(kind=8) :: get_DM_vectors(atom_all_shells,3,1)
! internal
integer :: i,natom,n_at_mag,n_at_nonmag,j,k,l
integer :: n_at_nonmag_incell,n_at_mag_inshell
real(kind=8) :: R1(3),R2(3),DMV(3),test,test_pos_A(3),test_pos_B(3),dumy
real(kind=8),allocatable :: R_mag(:,:),R_non_mag(:,:)
real(kind=8),allocatable :: pos_mag_atom_in_shell(:,:),pos_nonmag_atom_in_shell(:,:)
logical :: if_present,exists

natom=size(my_motif%atomic)
exists=.false.

n_at_mag=count(my_motif%i_mom(:),1)
n_at_nonmag=natom-n_at_mag
get_DM_vectors=0.0d0

allocate(R_mag(3,n_at_mag),R_non_mag(3,n_at_nonmag))
j=0
k=0
do i=1,natom
   if (abs(my_motif%atomic(i)%moment).gt.1.0d-8) then
      j=j+1
      R_mag(:,j)=matmul(transpose(areal),my_motif%atomic(i)%position(:))
   else
      k=k+1
      R_non_mag(:,k)=matmul(transpose(areal),my_motif%atomic(i)%position(:))
   endif

enddo

! find all magnetic atoms in the unit cell
! the number of directions for the different atoms are put in pos_mag_atom_in_shell
! for that apply symmetry operations
n_at_mag_inshell=find_all_equivalent_atom(R_mag(:,1)+areal(1,:),n_sym,symmetries,areal)
allocate(pos_mag_atom_in_shell(3,n_at_mag_inshell))
pos_mag_atom_in_shell=0.0d0
call find_all_atom(pos_mag_atom_in_shell,R_mag(:,1)+areal(1,:),n_sym,symmetries,areal)

! find all equivalent non magnetic atoms
! for that apply symmetry operations
n_at_nonmag_incell=find_all_equivalent_atom(R_non_mag(:,1)-R_mag(:,1),n_sym,symmetries,areal)
allocate(pos_nonmag_atom_in_shell(3,n_at_nonmag_incell))
call find_all_atom(pos_nonmag_atom_in_shell,R_non_mag(:,1)-R_mag(:,1),n_sym,symmetries,areal)

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

#ifdef CPP_DEBUG
do i=1,atom_all_shells
   write(*,*) get_DM_vectors(i,:,1)
enddo

#endif

end function get_DM_vectors












!
! function that finds ne number of equivalent atoms in the unit cell
!

function find_all_equivalent_atom(R,n_sym,symmetries,areal)
implicit none
real(kind=8), intent(in) :: R(:),areal(:,:)
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

subroutine find_all_atom(pos_atom_in_shell,R,n_sym,symmetries,areal)
implicit none
real(kind=8), intent(inout) :: pos_atom_in_shell(:,:)
real(kind=8), intent(in) :: R(:),areal(:,:)
integer, intent(in) :: n_sym
type(symop), intent(in) :: symmetries(:)
! internal
real(kind=8) :: test_pos(3)
integer :: j,l,k
logical :: if_present

j=2
pos_atom_in_shell(:,1)=R(:)

do k=1,n_sym

   test_pos=matmul(symmetries(k)%mat,R(:))

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
! routine that calculates the direction of the DM vector depending on the position of the non-magnetic atoms in the unit cell
! inputs
! ndm: number of DM vectors to calculate
! nei: number of shelves on which the DM are vectors are acting
! ind: index of the neighbors to consider
! r: lattice vectors
! motif: motif of atoms in the unit cell
! dim_lat: size of the supercell
! world: size of the world. The sze of the array matches the dimension of the problem: 1 for 1D, 2 for 2Ds...
! phase: 1 if thin films, 2 if multilayers
! c1: dummy to separate the cases
!
! output:
! a big array containing the DM vectors in cartesian coordinates


       function setup_DM_1D(ndm,r,motif,phase)
       use m_vector, only : cross,norm
       use m_sym_utils, only : order_zaxis,angle_oriented,rot_mat
       use m_constants, only : pi
       implicit none
! variable that come in
       integer, intent(in) :: ndm,phase
       real (kind=8), intent(in) :: r(3)
       type (cell), intent(in) :: motif
! value of the function
       real (kind=8) :: setup_DM_1D(ndm,3,phase)
! dummy variable
       integer :: i,j,i_dm,k
       real(kind=8) :: vec(3),R1(3),R2(3),DM_vec(3)
       real(kind=8) :: non_mag_at(count(.not.motif%i_mom),3)
! part of the symmetry

       i_DM=0
       setup_DM_1D=0
       non_mag_at=0.0d0

! find none magnetic atoms in the motif
       i=0
       do j=1,size(motif%i_mom)
        if (motif%i_mom(j)) cycle
        i=i+1
        non_mag_at(i,:)=motif%pos(j,:)
       enddo
      ! let's order them from the closest to the (0,0,0) atom to the furthest.
       do j=1,i
        do k=1,i
         if (abs(non_mag_at(j,3)).gt.abs(non_mag_at(k,3))) then
         vec=non_mag_at(j,:)
         non_mag_at(j,:)=non_mag_at(k,:)
         non_mag_at(k,:)=vec
         endif
        enddo
       enddo

       if (count(motif%i_mom).eq.1) then
        setup_DM_1D=0.0d0
       else
        do j=1,size(motif%i_mom)
         if (.not.motif%i_mom(j).and.(abs(sum(motif%pos(j,1:3))).lt.1.0d-8)) cycle
         do i=1,size(non_mag_at,1)
          R1=r(:)*non_mag_at(i,1)+non_mag_at(i,2)+non_mag_at(i,3)
          vec=r(:)*motif%pos(j,1)+motif%pos(j,2)+motif%pos(j,3)
          R2=R1+vec

          DM_vec=1/norm(R2)/norm(R1)/norm(R1-R2)*cross(R1,R2,1,3)
          setup_DM_1D(1,:,1)=DM_vec
          setup_DM_1D(2,2:3,1)=DM_vec(2:3)
          setup_DM_1D(2,1,1)=-DM_vec(1)
         enddo
        enddo
       endif

       end function setup_DM_1D

function setup_DM_2D(ndm,nei,ind,r,motif,world,phase)
use m_vector, only : cross,norm
use m_sym_utils, only : order_zaxis,angle_oriented,rot_mat,pos_nei
use m_constants, only : pi
use m_table_dist, only : Tdir
implicit none
! variable that come in
integer, intent(in) :: nei,ndm,world(:),phase
real (kind=8), intent(in) :: r(3,3)
integer, intent(in) :: ind(nei)
type (cell), intent(in) :: motif
! value of the function
real (kind=8) :: setup_DM_2D(ndm,3,phase)
! dummy variable
integer :: i,j,i_dm,n_z,k,step,i_nei,avant
real(kind=8) :: vec(3),R1(3),R1_dum(3),R2(3),dumy,DM_vec(3),vec_2(3)
real(kind=8) :: non_mag_at(count(.not.motif%i_mom),3),mag_at(count(motif%i_mom),3)
! directions of the neighbors
real(kind=8) :: tabledir(3,nei)
! part of the symmetry
real(kind=8) :: rot_ang,rotation(3,3)
logical :: exists

step=1
i_DM=0
setup_DM_2D=0
non_mag_at=0.0d0
mag_at=0.0d0
tabledir=0.0d0
avant=0

! find none magnetic atoms in the motif
i=0
k=0

do j=1,size(motif%i_mom)
   if (motif%i_mom(j))then
      k=k+1
      mag_at(k,:)=motif%atomic(j)%position
   else
      i=i+1
      non_mag_at(i,:)=motif%atomic(j)%position
   endif
enddo

call Tdir(r,nei,world,motif,tabledir)

n_z=order_zaxis(r)

        ! one magnetic atom in the unit cell, the rotation is done from u to v
if ((count(motif%i_mom).eq.1).or.(phase.ge.2)) then

   do i_nei=1,nei
      do i=1,size(non_mag_at,1)
         vec=r(1,:)*mag_at(1,1)+r(2,:)*mag_at(1,2)+r(3,:)*mag_at(1,3)
         vec_2=tabledir(:,i_nei)
         R1=r(1,:)*non_mag_at(i,1)+r(2,:)*non_mag_at(i,2)+r(3,:)*non_mag_at(i,3)
! find the position of the good non magnetic atom (goes to sym_utils.f90 for details)
         R1_dum=pos_nei(R1,vec,vec_2,r)

         R1=R1_dum-vec
         R2=R1_dum-vec_2
         DM_vec=cross(R1,R2,1,3)

         dumy=norm(DM_vec)
         DM_vec=DM_vec/dumy
             ! make sure that DM_vec and vec rotates in the sense of u toward v
!           write(*,*) 360.0d0/pi(2.0d0)*angle_oriented(vec,DM_vec)

           ! change by SvM: dble(j-1) -> dble((j-1)/2)
           ! The old version produced problems with the square lattice (0 degree, 180 degree, 360 degree=0 degree)
           ! while it was working for the hexagonal lattice (0 degree, 240 degree 480 degree=120 degree)
         do j=1,ind(i_nei),2
            rot_ang=360.0d0/dble(n_z)*dble((j-1)/2)
            rotation=rot_mat(rot_ang,(/0,0,1/))
            setup_DM_2D(avant+j,:,i)=matmul(rotation,DM_vec)
         enddo

!           do j=1,ind(i_nei),2
!            rot_ang=360.0d0/dble(n_z)*dble(j-1)
!            rotation=rot_mat(rot_ang,(/0,0,1/))
!            setup_DM_2D(avant+j,:,i)=matmul(rotation,DM_vec)
!           enddo

         do j=2,ind(i_nei),2
            setup_DM_2D(avant+j,:,i)=-setup_DM_2D(avant+j-1,:,i)
         enddo

      enddo
      avant=avant+ind(i_nei)
   enddo

!!!!!!!!!!
!!!!!!!!!!
elseif (count(motif%i_mom).ne.1) then
    write(6,'(a)') "not coded"
endif

do j=1,phase
   do i=1,ndm
      dumy=norm(setup_DM_2D(i,:,j))
      if (dumy.gt.1.0d-8) setup_DM_2D(i,:,j)=setup_DM_2D(i,:,j)/dumy
   enddo
enddo

inquire(file='DM-2donly',exist=exists)
if (exists) then
   setup_DM_2D(:,3,:)=0.0d0
   do j=1,phase
      do i=1,ndm
         dumy=norm(setup_DM_2D(i,:,j))
         setup_DM_2D(i,:,j)=setup_DM_2D(i,:,j)/dumy
      enddo
   enddo
endif

#ifdef CPP_DEBUG
       write(*,*) ndm,nei,phase
       do j=1,phase
       i=1
       do while (i.le.ndm)
       write(*,*) setup_DM_2D(i,:,j)
       i=i+1
       enddo
       enddo
#endif

end function setup_DM_2D

       function setup_DM_3D(ndm,r,motif,phase)
       use m_vector, only : cross,norm
       use m_sym_utils, only : order_zaxis,angle_oriented,rot_mat
       use m_constants, only : pi
       implicit none
! variable that come in
       integer, intent(in) :: ndm,phase
       real (kind=8), intent(in) :: r(3,3)
       type (cell), intent(in) :: motif
! value of the function
       real (kind=8) :: setup_DM_3D(ndm,3,phase)
! dummy variable
       integer :: i,j,i_dm,k
       real(kind=8) :: vec(3),R1(3),R2(3),DM_vec(3)
       real(kind=8) :: non_mag_at(count(.not.motif%i_mom),3)

       i_DM=0
       setup_DM_3D=0
       non_mag_at=0.0d0

! find none magnetic atoms in the motif
       i=0
       do j=1,size(motif%i_mom)
        if (motif%i_mom(j)) cycle
        i=i+1
        non_mag_at(i,:)=motif%pos(j,:)
       enddo
      ! let's order them from the closest to the (0,0,0) atom to the furthest.
       do j=1,i
        do k=1,i
         if (abs(non_mag_at(j,3)).gt.abs(non_mag_at(k,3))) then
         vec=non_mag_at(j,:)
         non_mag_at(j,:)=non_mag_at(k,:)
         non_mag_at(k,:)=vec
         endif
        enddo
       enddo

       if (count(motif%i_mom).eq.1) then
        setup_DM_3D=0.0d0
       else
        do j=1,size(motif%i_mom)
         if (.not.motif%i_mom(j).and.(abs(sum(motif%pos(j,1:3))).lt.1.0d-8)) cycle
         do i=1,size(non_mag_at,1)
          R1=r(1,:)*non_mag_at(i,1)+r(2,:)*non_mag_at(i,2)+r(3,:)*non_mag_at(i,3)
          vec=r(1,:)*motif%pos(j,1)+r(2,:)*motif%pos(j,2)+r(3,:)*motif%pos(j,3)
          R2=R1+vec

          DM_vec=1/norm(R2)/norm(R1)/norm(R1-R2)*cross(R1,R2,1,3)
          setup_DM_3D(1,:,1)=DM_vec
          setup_DM_3D(2,2:3,1)=DM_vec(2:3)
          setup_DM_3D(2,1,1)=-DM_vec(1)
         enddo
        enddo
       endif

       end function setup_DM_3D

       end module
