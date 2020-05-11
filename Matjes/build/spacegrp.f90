!#######  COPYRIGHT ########
!#  Copyright 1995-2014 Timo Schena, Phivos Mavropoulos, Yurij Mokrousov, Cyrille Barreteau and Stefan Bluegel
!#
!#    This file is part of JuTiBi.
!#    JuTiBi is free software: you can redistribute it and/or modify
!#    it under the terms of the GNU General Public License as published by
!#    the Free Software Foundation, either version 3 of the License, or
!#    (at your option) any later version.
!#    JuTiBi is distributed in the hope that it will be useful,
!#    but WITHOUT ANY WARRANTY; without even the implied warranty of
!#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!#    GNU General Public License for more details.
!#    You should have received a copy of the GNU General Public License
!#    along with JuTiBi.  If not, see <http://www.gnu.org/licenses/>.
!###########################

module m_spacegrp
use m_derived_types

private
public :: findgroup
contains

subroutine findgroup(bravais,bravais1,recbv,recbv1,
     &                     rbasis,rbasis1,alat, abclat ,nbasis,
     &                     rsymat,rotname,nsymat,isymindex)

implicit none
! **********************************************************
! This subroutine finds the rotation matrices that leave the
! real lattice unchanged.
! input:  bravais(i,j)    true bravais lattice vectors
!                         i = x,y,z ; j = A, B, C in units of (a,b,c)
!         recbv(i,j)      reciprocal basis vectors in units of (1/a,1/b,1/c)
!         rbasis1         coordinates of basis atoms in units of (a,b,c)
!         nbasis          number of basis atoms
!         alat            lattice constants
!         rsymat          all 64 rotation matrices.
!         rotname         names for the rotation matrices
! output: nsymat          number of rotations that restore the lattice.
!         ISYMINDEX       index for the symmetries found
!         bravais1        bravais vectors in angstrom (or a.u.)
!         recbv1          reciprocal bravais vectors in 1/angstrom (or 1/a.u.)
!         rbasis
!
! This sub makes all 64 rotations in the basis vectors and bravais
! vectors and checks if the new rotated vector belongs in the
! lattice. The proper rotation must bring all vectors to a lattice
! vector. Information about the rotations found is printed in the end.
! The array ISYMINDEX holds the numbers of the symmetry operations
! that are stored in array RSYMAT
! **********************************************************

integer :: NSYMAXD
integer, parameter :: NSYMAXD=48
integer :: nbasis,nsymat
integer :: isymindex(NSYMAXD)
real(kind=8) :: BRAVAIS(3,3)
real(kind=8), dimension(3,nbasis) :: RBASIS , RBASIS1
real(kind=8) RSYMAT(64,3,3),recbv(3,3),recbv1(3,3)
!
! Local variables
!
real(kind=8) :: r(3,4)
real(kind=8), dimension(3,nbasis) :: rotrbas
real(kind=8) ::  alat(3), abclat(3) ,bravais1(3,3) !abclat contains a, b , c as lattice constants not a , b/a and c/a
integer :: i,j,isym,nsym,ia
logical llatbas,LBULK

!     -------------------------------------------------------------
      NSYM = 0
!     - ---------------------------------

write(6,'(a)') '################################'
write(6,'(a)') 'findgroup: Calculating Bravais-',
     &           ' and Basisvectors into a.u.'

      abclat(1)=alat(1)
      abclat(2)=alat(2)*alat(1)
      abclat(3)=alat(3)*alat(1)

      do i=1,3
         do j=1,3
            bravais1(j,i) = bravais(j,i)*abclat(j)
            recbv1(j,i) = recbv(j,i)/abclat(j)
         end do
      end do


      do i=1,nbasis
         do j=1,3
            RBASIS1(j,i)=RBASIS(j,i)*abclat(j)
         end do
      end do

!     Check for surface mode. If so, set bravais1(3,3) very large, so
!     that only the in-plane symmetries are found. Not checked, be careful of z--> -z!
      LBULK=.TRUE.
!     Now check the bravais vectors if they have a z component
      if ((bravais(1,3).eq.0.d0).and.(bravais(2,3).eq.0.d0).and.
     &     (bravais(3,3).eq.0.d0)) THEN
         LBULK=.FALSE.
      END IF
!
      do isym=1,64
!     rotate bravais lattice vectors
!
!     In the case of slab/interface geometry look only for
!     symmetry opperations that preserve the z axis..
!
         IF (LBULK .OR. (RSYMAT(ISYM,3,3).EQ.1) ) THEN
!     do rotation only in case bulk or if slab and z axis is restored..


            do i=1,3            ! Loop on bravais vectors
               do j=1,3         ! Loop on coordinates
                  r(j,i) = rsymat(isym,j,1)*bravais1(1,i) +
     &                 rsymat(isym,j,2)*bravais1(2,i) +
     &                 rsymat(isym,j,3)*bravais1(3,i)
               enddo
            enddo
!
!     rotate the basis atoms p and take RSYMAT.p - p then
!     find if R = (RSYMAT.bravais + RSYMAT.p - p) belongs to the
!     lattice. This is done by function latvec by checking
!     if R.q = integer (q reciprocal lattice vector)
!
            llatbas = .true.
            do ia=1,nbasis      ! Loop on basis atoms
               do j=1,3         ! Loop on coordinates
                  rotrbas(j,ia) = rsymat(isym,j,1)*rbasis1(1,ia) +
     &                 rsymat(isym,j,2)*rbasis1(2,ia) +
     &                 rsymat(isym,j,3)*rbasis1(3,ia)
!
                  rotrbas(j,ia) = rotrbas(j,ia) - rbasis1(j,ia)
                  r(j,4) = rotrbas(j,ia)   ! here are the basis vectors for the symmetry examination
               enddo
               if (.not.latvec(4,recbv1,r)) llatbas=.false. !attention: basis atoms aren't treat equivalent!!
            enddo               ! ia=1,nbasis

!
!     if llatbas=.true. the rotation does not change the lattice
!
            if (llatbas) then
               NSYM = NSYM + 1
               ISYMINDEX(NSYM) = ISYM
            end if
         END IF                 ! (LBULK .OR. (RSYMAT(ISYM,3,3).EQ.1) )
      end do                    ! isym=1,nmatd
!     nsym symmetries were found
!     the ISYMINDEX array has the numbers of the symmetries found
!
!
      NSYMAT = NSYM
      WRITE(*,*) 'findgroup: Number of symmetries: ', NSYMAT
      write(*,*) 'findgroup: Symmetry-operations: ',
     &             (rotname(isymindex(isym)),isym=1,NSYMAT)

end subroutine findgroup

logical function latvec(n,qlat,vec)
!- Checks if a set of vectors are lattice vectors
! ----------------------------------------------------------------------
!i Inputs:
!i   n     :number of vectors
!i   qlat  :primitive translation vectors in reciprocal space
!i   vec   :double-precision vector
!o Outputs:
!o   latvec:.true. if all vectors are lattice vectors
!r Remarks:
! ----------------------------------------------------------------------
implicit none
! Passed parameters:
integer :: n
real(kind=8) :: qlat(3,*),vec(3,*)
! Local parameters:
integer :: i,m
real(kind=8) :: vdiff
real(kind=8),parameter :: tol=1.d-6
! Common block:
! Intrinsic functions:
intrinsic  dabs,dnint

latvec=.false.
do i=1,n
   do m=1,3
      vdiff=vec(1,i)*qlat(1,m)+vec(2,i)*qlat(2,m)+vec(3,i)*qlat(3,m)
      vdiff=dabs(vdiff-dnint(vdiff))
      if (vdiff.gt.tol) return
   enddo
enddo
latvec=.true.
END function latvec

end module m_spacegrp

