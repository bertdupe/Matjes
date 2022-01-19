module m_symmetry_base
use m_basic_types, only : symop

implicit none

private
public : pt_grp


! type symmetry operation
type, abstract :: pt_grp
 type(symop),allocatable :: rotmat(:)
 integer                 :: n_sym=0
end type

 contains
 procedure(init_init),deferred :: init
 procedure(apply)    ,deferred :: apply_sym

 procedure,NON_OVERRIDABLE   :: init_base
 procedure,NON_OVERRIDABLE   :: load_base
end type symop

abstract interface
   subroutine int_init(this,N)
        import pt_grp
        class(pt_grp),intent(out)       :: this
        integer      ,intent(in)        :: N_sym
    end subroutine

   function apply(this,i,u) result(v)
        import pt_grp
        class(pt_grp),intent(in)       :: this
        integer      ,intent(in)       :: i
        real(8)      ,intent(in)       :: u(3)
        real(8)                        :: v(3)
    end function
end interface






contains




subroutine init_base(this)
type(this),intent(out) :: this(:)
integer i

allocate(this(64))
! set to zero
do i=1,64
     this(i1)%mat = 0.d0
     this(i1)%name=""
end do
this%n_sym=64

end subroutine

subroutine load_base(this)
! **********************************************
! This subroutine defines the rotation matrices for
! all the 32 point groups and names them after
! J.F. Cornwell (Group Theory??) second edition
! Appendix D, p 324-325
!
!  rotmat , real*8(3,3,64) : All 64 3x3 rotation matriced of the point groups
!  rotname , char*10(64)   : Names of the 64 rotation matrixes after J.F. Cornwell
!
! the basis vectors have the fomat (line,column)
! *********************************************
implicit none
type(this),intent(out) :: this(:)
integer i1,is
real(kind=8) :: RTHREE,HALF

RTHREE = SQRT(3.d0)/2.d0
HALF = 0.5d0

!
      this(1)%mat(1,1) =  1.d0
      this(1)%mat(2,2) =  1.d0
      this(1)%mat(3,3) =  1.d0
      this(1)%name = 'E'
!
      this(2)%mat(1,2) =  1.d0
      this(2)%mat(2,3) = -1.d0
      this(2)%mat(3,1) = -1.d0
      this(2)%name = 'C3alfa'
!
      this(3)%mat(1,2) = -1.d0
      this(3)%mat(2,3) = -1.d0
      this(3)%mat(3,1) =  1.d0
      this(3)%name = 'C3beta '
!
      this(4)%mat(1,2) = -1.d0
      this(4)%mat(2,3) =  1.d0
      this(4)%mat(3,1) = -1.d0
      this(4)%name = 'C3gamma'
!
      this(5)%mat(1,2) = 1.d0
      this(5)%mat(2,3) = 1.d0
      this(5)%mat(3,1) = 1.d0
      this(5)%name = 'C3delta '
!
      this(6)%mat(1,3) = -1.d0
      this(6)%mat(2,1) =  1.d0
      this(6)%mat(3,2) = -1.d0
      this(6)%name = 'C3alfa-1'
!
      this(7)%mat(1,3) =  1.d0
      this(7)%mat(2,1) = -1.d0
      this(7)%mat(3,2) = -1.d0
      this(7)%name = 'C3beta-1 '
!
      this(8)%mat(1,3) = -1.d0
      this(8)%mat(2,1) = -1.d0
      this(8)%mat(3,2) =  1.d0
      this(8)%name = 'C3gamma-1'
!
      this(9)%mat(1,3) =  1.d0
      this(9)%mat(2,1) =  1.d0
      this(9)%mat(3,2) =  1.d0
      this(9)%name = 'C3delta-1'
!
      this(10)%mat(1,1) =  1.d0
      this(10)%mat(2,2) = -1.d0
      this(10)%mat(3,3) = -1.d0
      this(10)%name = 'C2x'
!
      this(11)%mat(1,1) = -1.d0
      this(11)%mat(2,2) =  1.d0
      this(11)%mat(3,3) = -1.d0
      this(11)%name = 'C2y'
!
      this(12)%mat(1,1) = -1.d0
      this(12)%mat(2,2) = -1.d0
      this(12)%mat(3,3) =  1.d0
      this(12)%name = 'C2z'
!
      this(13)%mat(1,1) =  1.d0
      this(13)%mat(2,3) =  1.d0
      this(13)%mat(3,2) = -1.d0
      this(13)%name = 'C4x'
!
      this(14)%mat(1,3) = -1.d0
      this(14)%mat(2,2) =  1.d0
      this(14)%mat(3,1) =  1.d0
      this(14)%name = 'C4y '
!
      this(15)%mat(1,2) =  1.d0
      this(15)%mat(2,1) = -1.d0
      this(15)%mat(3,3) =  1.d0
      this(15)%name = 'C4z'
!
      this(16)%mat(1,1) =  1.d0
      this(16)%mat(2,3) = -1.d0
      this(16)%mat(3,2) =  1.d0
      this(16)%name = 'C4x-1 '
!
      this(17)%mat(1,3) =  1.d0
      this(17)%mat(2,2) =  1.d0
      this(17)%mat(3,1) = -1.d0
      this(17)%name = 'C4y-1'
!
      this(18)%mat(1,2) = -1.d0
      this(18)%mat(2,1) =  1.d0
      this(18)%mat(3,3) =  1.d0
      this(18)%name = 'C4z-1'
!
      this(19)%mat(1,2) =  1.d0
      this(19)%mat(2,1) =  1.d0
      this(19)%mat(3,3) = -1.d0
      this(19)%name = 'C2a'
!
      this(20)%mat(1,2) = -1.d0
      this(20)%mat(2,1) = -1.d0
      this(20)%mat(3,3) = -1.d0
      this(20)%name = 'C2b'
!
      this(21)%mat(1,3) =  1.d0
      this(21)%mat(2,2) = -1.d0
      this(21)%mat(3,1) =  1.d0
      this(21)%name = 'C2c'
!
      this(22)%mat(1,3) = -1.d0
      this(22)%mat(2,2) = -1.d0
      this(22)%mat(3,1) = -1.d0
      this(22)%name = 'C2d'
!
      this(23)%mat(1,1) = -1.d0
      this(23)%mat(2,3) =  1.d0
      this(23)%mat(3,2) =  1.d0
      this(23)%name = 'C2e'
!
      this(24)%mat(1,1) = -1.d0
      this(24)%mat(2,3) = -1.d0
      this(24)%mat(3,2) = -1.d0
      this(24)%name = 'C2f'

do i1=1,24
   this(i1+24)%mat = -this(i1)%mat
   this(i1+24)%name = 'I'//this(i1)%name(1:9)
end do
!
!
!*********************************************
! Trigonal and hexagonal groups
!*********************************************
!
      this(49)%mat(1,1) = -HALF
      this(49)%mat(1,2) =  RTHREE
      this(49)%mat(2,1) = -RTHREE
      this(49)%mat(2,2) = -HALF
      this(49)%mat(3,3) =  1.d0
      this(49)%name = 'C3z'
!
      this(50)%mat(1,1) = -HALF
      this(50)%mat(1,2) = -RTHREE
      this(50)%mat(2,1) =  RTHREE
      this(50)%mat(2,2) = -HALF
      this(50)%mat(3,3) =  1.d0
      this(50)%name = 'C3z-1'
!
      this(51)%mat(1,1) =  HALF
      this(51)%mat(1,2) =  RTHREE
      this(51)%mat(2,1) = -RTHREE
      this(51)%mat(2,2) =  HALF
      this(51)%mat(3,3) =  1.d0
      this(51)%name = 'C6z'
!
      this(52)%mat(1,1) =  HALF
      this(52)%mat(1,2) = -RTHREE
      this(52)%mat(2,1) =  RTHREE
      this(52)%mat(2,2) =  HALF
      this(52)%mat(3,3) =  1.d0
      this(52)%name = 'C6z-1'
!
      this(53)%mat(1,1) = -HALF
      this(53)%mat(1,2) =  RTHREE
      this(53)%mat(2,1) =  RTHREE
      this(53)%mat(2,2) =  HALF
      this(53)%mat(3,3) = -1.d0
      this(53)%name = 'C2A'
!
      this(54)%mat(1,1) = -HALF
      this(54)%mat(1,2) = -RTHREE
      this(54)%mat(2,1) = -RTHREE
      this(54)%mat(2,2) =  HALF
      this(54)%mat(3,3) = -1.d0
      this(54)%name = 'C2B'
!
      this(55)%mat(1,1) =  HALF
      this(55)%mat(1,2) = -RTHREE
      this(55)%mat(2,1) = -RTHREE
      this(55)%mat(2,2) = -HALF
      this(55)%mat(3,3) = -1.d0
      this(55)%name = 'C2C'
!
      this(56)%mat(1,1) =  HALF
      this(56)%mat(1,2) =  RTHREE
      this(56)%mat(2,1) =  RTHREE
      this(56)%mat(2,2) = -HALF
      this(56)%mat(3,3) = -1.d0
      this(56)%name = 'C2D'

do is=1,8
   this(56+is)%mat = -ROTMAT(48+is)%mat
   this(56+is)%name = 'I'//ROTMAT(48+is)%name(1:9)
end do

!-----------------------------------------
end subroutine

end module

end module
