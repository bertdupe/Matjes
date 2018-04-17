module m_pointgrp
use m_basic_types, only : symop

private
public :: get_symop
contains

subroutine get_symop(rotmat)
! **********************************************
! This subroutine defines the rotation matrices for
! all the 32 point groups and names them after
! J.F. Cornwell (Group Theory??) second edition
! Appendix D, p 324-325
!
!  rotmat , real*8(3,3,64) : All 64 3x3 rotation matriced of the point groups
!  rotname , char*10(64)   : Names of the 64 rotation matrixes after J.F. Cornwell
!
! *********************************************
implicit none
integer i1,is
type(symop),intent(out) :: ROTMAT(64)
real(kind=8) :: RTHREE,HALF

RTHREE = SQRT(3.d0)/2.d0
HALF = 0.5d0
! set to zero
do i1=1,64
     ROTMAT(i1)%mat = 0.d0
end do
!
      ROTMAT(1)%mat(1,1) =  1.d0
      ROTMAT(1)%mat(2,2) =  1.d0
      ROTMAT(1)%mat(3,3) =  1.d0
      ROTMAT(1)%name = 'E'
!
      ROTMAT(2)%mat(1,2) =  1.d0
      ROTMAT(2)%mat(2,3) = -1.d0
      ROTMAT(2)%mat(3,1) = -1.d0
      ROTMAT(2)%name = 'C3alfa'
!
      ROTMAT(3)%mat(1,2) = -1.d0
      ROTMAT(3)%mat(2,3) = -1.d0
      ROTMAT(3)%mat(3,1) =  1.d0
      ROTMAT(3)%name = 'C3beta '
!
      ROTMAT(4)%mat(1,2) = -1.d0
      ROTMAT(4)%mat(2,3) =  1.d0
      ROTMAT(4)%mat(3,1) = -1.d0
      ROTMAT(4)%name = 'C3gamma'
!
      ROTMAT(5)%mat(1,2) = 1.d0
      ROTMAT(5)%mat(2,3) = 1.d0
      ROTMAT(5)%mat(3,1) = 1.d0
      ROTMAT(5)%name = 'C3delta '
!
      ROTMAT(6)%mat(1,3) = -1.d0
      ROTMAT(6)%mat(2,1) =  1.d0
      ROTMAT(6)%mat(3,2) = -1.d0
      ROTMAT(6)%name = 'C3alfa-1'
!
      ROTMAT(7)%mat(1,3) =  1.d0
      ROTMAT(7)%mat(2,1) = -1.d0
      ROTMAT(7)%mat(3,2) = -1.d0
      ROTMAT(7)%name = 'C3beta-1 '
!
      ROTMAT(8)%mat(1,3) = -1.d0
      ROTMAT(8)%mat(2,1) = -1.d0
      ROTMAT(8)%mat(3,2) =  1.d0
      ROTMAT(8)%name = 'C3gamma-1'
!
      ROTMAT(9)%mat(1,3) =  1.d0
      ROTMAT(9)%mat(2,1) =  1.d0
      ROTMAT(9)%mat(3,2) =  1.d0
      ROTMAT(9)%name = 'C3delta-1'
!
      ROTMAT(10)%mat(1,1) =  1.d0
      ROTMAT(10)%mat(2,2) = -1.d0
      ROTMAT(10)%mat(3,3) = -1.d0
      ROTMAT(10)%name = 'C2x'
!
      ROTMAT(11)%mat(1,1) = -1.d0
      ROTMAT(11)%mat(2,2) =  1.d0
      ROTMAT(11)%mat(3,3) = -1.d0
      ROTMAT(11)%name = 'C2y'
!
      ROTMAT(12)%mat(1,1) = -1.d0
      ROTMAT(12)%mat(2,2) = -1.d0
      ROTMAT(12)%mat(3,3) =  1.d0
      ROTMAT(12)%name = 'C2z'
!
      ROTMAT(13)%mat(1,1) =  1.d0
      ROTMAT(13)%mat(2,3) =  1.d0
      ROTMAT(13)%mat(3,2) = -1.d0
      ROTMAT(13)%name = 'C4x'
!
      ROTMAT(14)%mat(1,3) = -1.d0
      ROTMAT(14)%mat(2,2) =  1.d0
      ROTMAT(14)%mat(3,1) =  1.d0
      ROTMAT(14)%name = 'C4y '
!
      ROTMAT(15)%mat(1,2) =  1.d0
      ROTMAT(15)%mat(2,1) = -1.d0
      ROTMAT(15)%mat(3,3) =  1.d0
      ROTMAT(15)%name = 'C4z'
!
      ROTMAT(16)%mat(1,1) =  1.d0
      ROTMAT(16)%mat(2,3) = -1.d0
      ROTMAT(16)%mat(3,2) =  1.d0
      ROTMAT(16)%name = 'C4x-1 '
!
      ROTMAT(17)%mat(1,3) =  1.d0
      ROTMAT(17)%mat(2,2) =  1.d0
      ROTMAT(17)%mat(3,1) = -1.d0
      ROTMAT(17)%name = 'C4y-1'
!
      ROTMAT(18)%mat(1,2) = -1.d0
      ROTMAT(18)%mat(2,1) =  1.d0
      ROTMAT(18)%mat(3,3) =  1.d0
      ROTMAT(18)%name = 'C4z-1'
!
      ROTMAT(19)%mat(1,2) =  1.d0
      ROTMAT(19)%mat(2,1) =  1.d0
      ROTMAT(19)%mat(3,3) = -1.d0
      ROTMAT(19)%name = 'C2a'
!
      ROTMAT(20)%mat(1,2) = -1.d0
      ROTMAT(20)%mat(2,1) = -1.d0
      ROTMAT(20)%mat(3,3) = -1.d0
      ROTMAT(20)%name = 'C2b'
!
      ROTMAT(21)%mat(1,3) =  1.d0
      ROTMAT(21)%mat(2,2) = -1.d0
      ROTMAT(21)%mat(3,1) =  1.d0
      ROTMAT(21)%name = 'C2c'
!
      ROTMAT(22)%mat(1,3) = -1.d0
      ROTMAT(22)%mat(2,2) = -1.d0
      ROTMAT(22)%mat(3,1) = -1.d0
      ROTMAT(22)%name = 'C2d'
!
      ROTMAT(23)%mat(1,1) = -1.d0
      ROTMAT(23)%mat(2,3) =  1.d0
      ROTMAT(23)%mat(3,2) =  1.d0
      ROTMAT(23)%name = 'C2e'
!
      ROTMAT(24)%mat(1,1) = -1.d0
      ROTMAT(24)%mat(2,3) = -1.d0
      ROTMAT(24)%mat(3,2) = -1.d0
      ROTMAT(24)%name = 'C2f'

do i1=1,24
   ROTMAT(i1+24)%mat = -ROTMAT(i1)%mat
   ROTMAT(i1+24)%name = 'I'//ROTMAT(i1)%name(1:9)
end do
!
!
!*********************************************
! Trigonal and hexagonal groups
!*********************************************
!
      ROTMAT(49)%mat(1,1) = -HALF
      ROTMAT(49)%mat(1,2) =  RTHREE
      ROTMAT(49)%mat(2,1) = -RTHREE
      ROTMAT(49)%mat(2,2) = -HALF
      ROTMAT(49)%mat(3,3) =  1.d0
      ROTMAT(49)%name = 'C3z'
!
      ROTMAT(50)%mat(1,1) = -HALF
      ROTMAT(50)%mat(1,2) = -RTHREE
      ROTMAT(50)%mat(2,1) =  RTHREE
      ROTMAT(50)%mat(2,2) = -HALF
      ROTMAT(50)%mat(3,3) =  1.d0
      ROTMAT(50)%name = 'C3z-1'
!
      ROTMAT(51)%mat(1,1) =  HALF
      ROTMAT(51)%mat(1,2) =  RTHREE
      ROTMAT(51)%mat(2,1) = -RTHREE
      ROTMAT(51)%mat(2,2) =  HALF
      ROTMAT(51)%mat(3,3) =  1.d0
      ROTMAT(51)%name = 'C6z'
!
      ROTMAT(52)%mat(1,1) =  HALF
      ROTMAT(52)%mat(1,2) = -RTHREE
      ROTMAT(52)%mat(2,1) =  RTHREE
      ROTMAT(52)%mat(2,2) =  HALF
      ROTMAT(52)%mat(3,3) =  1.d0
      ROTMAT(52)%name = 'C6z-1'
!
      ROTMAT(53)%mat(1,1) = -HALF
      ROTMAT(53)%mat(1,2) =  RTHREE
      ROTMAT(53)%mat(2,1) =  RTHREE
      ROTMAT(53)%mat(2,2) =  HALF
      ROTMAT(53)%mat(3,3) = -1.d0
      ROTMAT(53)%name = 'C2A'
!
      ROTMAT(54)%mat(1,1) = -HALF
      ROTMAT(54)%mat(1,2) = -RTHREE
      ROTMAT(54)%mat(2,1) = -RTHREE
      ROTMAT(54)%mat(2,2) =  HALF
      ROTMAT(54)%mat(3,3) = -1.d0
      ROTMAT(54)%name = 'C2B'
!
      ROTMAT(55)%mat(1,1) =  HALF
      ROTMAT(55)%mat(1,2) = -RTHREE
      ROTMAT(55)%mat(2,1) = -RTHREE
      ROTMAT(55)%mat(2,2) = -HALF
      ROTMAT(55)%mat(3,3) = -1.d0
      ROTMAT(55)%name = 'C2C'
!
      ROTMAT(56)%mat(1,1) =  HALF
      ROTMAT(56)%mat(1,2) =  RTHREE
      ROTMAT(56)%mat(2,1) =  RTHREE
      ROTMAT(56)%mat(2,2) = -HALF
      ROTMAT(56)%mat(3,3) = -1.d0
      ROTMAT(56)%name = 'C2D'

do is=1,8
   ROTMAT(56+is)%mat = -ROTMAT(48+is)%mat
   ROTMAT(56+is)%name = 'I'//ROTMAT(48+is)%name(1:9)
end do
      
!-----------------------------------------
end subroutine

end module m_pointgrp
