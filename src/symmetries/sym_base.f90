module m_symmetry_base
use m_basic_types, only : symop

implicit none

private
public :: pt_grp


! type symmetry operation
type,extends(symop), abstract :: pt_grp
   type(symop),allocatable :: rotmat(:)
   integer                 :: n_sym=0
   real(8)                 :: tol_sym=1.0d-6

 contains
   procedure(init_init),deferred :: init
   procedure(init_load),deferred :: load
   procedure(apply)    ,deferred :: apply_sym
   procedure(num_sym)  ,deferred :: get_N_sym
   procedure(latt_sym) ,deferred :: get_latt_sym
   procedure(pt_sym)   ,deferred :: get_pt_sym
   procedure(rw_sym)   ,deferred :: write_sym
   procedure(rd_sym)   ,deferred :: read_sym
   procedure(load_all) ,deferred :: get_all_symetries

   procedure,NON_OVERRIDABLE   :: init_base
   procedure,NON_OVERRIDABLE   :: load_base
   procedure,NON_OVERRIDABLE   :: read_test
end type



abstract interface
   subroutine init_init(this,N_sym)
        import pt_grp
        class(pt_grp),intent(out)       :: this
        integer      ,intent(in)        :: N_sym
    end subroutine

   subroutine init_load(this,index_sym,all_sym,sym_translation)
        import pt_grp
        class(pt_grp),intent(inout)     :: this
        class(pt_grp),intent(in)        :: all_sym
        integer      ,intent(in)        :: index_sym(:)
        real(8)      ,intent(in)        :: sym_translation(:,:)
    end subroutine

   function apply(this,i,u) result(v)
        import pt_grp
        class(pt_grp),intent(in)       :: this
        integer      ,intent(in)       :: i
        real(8)      ,intent(in)       :: u(3)
        real(8)                        :: v(3)
    end function

    function num_sym(this) result(N)
        import pt_grp
        class(pt_grp),intent(in)       :: this
        integer                        :: N
    end function

    subroutine latt_sym(this,areal,number_sym,sym_index,periodic,dim_lat)
        import pt_grp
        class(pt_grp),intent(in)       :: this
        real(8)      ,intent(in)       :: areal(3,3)
        integer      ,intent(out)      :: number_sym
        integer      ,intent(inout)    :: sym_index(:)
        integer      ,intent(in)       :: dim_lat(:)
        logical      ,intent(in)       :: periodic(:)
    end subroutine

    subroutine pt_sym(this,my_lattice,number_sym,sym_index,my_motif,sym_translation)
        use m_derived_types, only : t_cell
        use m_type_lattice , only : lattice
        import pt_grp
        class(pt_grp),intent(in)          :: this
        type(lattice), intent(in)         :: my_lattice
        integer      , intent(inout)      :: number_sym
        integer,intent(inout),allocatable :: sym_index(:)
        type(t_cell) ,intent(in)          :: my_motif
        real(8),intent(inout),allocatable :: sym_translation(:,:)
    end subroutine

    subroutine rw_sym(this)
        import pt_grp
        class(pt_grp),intent(in)       :: this
    end subroutine

    subroutine rd_sym(this,fname)
        import pt_grp
        class(pt_grp),intent(inout)    :: this
        character(len=*),intent(in)    :: fname
    end subroutine

!    subroutine rd_test(this,test)
!        import pt_grp
!        class(pt_grp),intent(inout)    :: this
!        logical      ,intent(out)      :: test
!    end subroutine

    subroutine load_all(this,N,out)
        import pt_grp
        class(pt_grp),intent(in)       :: this
        integer      ,intent(out)      :: N
        class(pt_grp),intent(inout)    :: out
    end subroutine

end interface






contains


subroutine read_test(this,test)
    use m_io_files_utils
    class(pt_grp)     , intent(inout) :: this
    logical           , intent(out)   :: test

    inquire(file='symmetries.in',exist=test)
    if (test) then
       write(6,'(a)') 'getting symmetries from file'
       call this%read_sym('symmetries.in')
       call this%write_sym()
    else
       write(6,'(a)') 'symmetries.in not found - calculating symmetries'
    endif

end subroutine

subroutine init_base(this)
class(pt_grp),intent(out) :: this
integer i

allocate(this%rotmat(64))
! set to zero
do i=1,64
     this%rotmat(i)%mat = 0.d0
     this%rotmat(i)%name=""
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
class(pt_grp),intent(inout) :: this
integer i1,is
real(kind=8) :: RTHREE,HALF

RTHREE = SQRT(3.d0)/2.d0
HALF = 0.5d0

!
      this%rotmat(1)%mat(1,1) =  1.d0
      this%rotmat(1)%mat(2,2) =  1.d0
      this%rotmat(1)%mat(3,3) =  1.d0
      this%rotmat(1)%name = 'E'
!
      this%rotmat(2)%mat(1,2) =  1.d0
      this%rotmat(2)%mat(2,3) = -1.d0
      this%rotmat(2)%mat(3,1) = -1.d0
      this%rotmat(2)%name = 'C3alfa'
!
      this%rotmat(3)%mat(1,2) = -1.d0
      this%rotmat(3)%mat(2,3) = -1.d0
      this%rotmat(3)%mat(3,1) =  1.d0
      this%rotmat(3)%name = 'C3beta '
!
      this%rotmat(4)%mat(1,2) = -1.d0
      this%rotmat(4)%mat(2,3) =  1.d0
      this%rotmat(4)%mat(3,1) = -1.d0
      this%rotmat(4)%name = 'C3gamma'
!
      this%rotmat(5)%mat(1,2) = 1.d0
      this%rotmat(5)%mat(2,3) = 1.d0
      this%rotmat(5)%mat(3,1) = 1.d0
      this%rotmat(5)%name = 'C3delta '
!
      this%rotmat(6)%mat(1,3) = -1.d0
      this%rotmat(6)%mat(2,1) =  1.d0
      this%rotmat(6)%mat(3,2) = -1.d0
      this%rotmat(6)%name = 'C3alfa-1'
!
      this%rotmat(7)%mat(1,3) =  1.d0
      this%rotmat(7)%mat(2,1) = -1.d0
      this%rotmat(7)%mat(3,2) = -1.d0
      this%rotmat(7)%name = 'C3beta-1 '
!
      this%rotmat(8)%mat(1,3) = -1.d0
      this%rotmat(8)%mat(2,1) = -1.d0
      this%rotmat(8)%mat(3,2) =  1.d0
      this%rotmat(8)%name = 'C3gamma-1'
!
      this%rotmat(9)%mat(1,3) =  1.d0
      this%rotmat(9)%mat(2,1) =  1.d0
      this%rotmat(9)%mat(3,2) =  1.d0
      this%rotmat(9)%name = 'C3delta-1'
!
      this%rotmat(10)%mat(1,1) =  1.d0
      this%rotmat(10)%mat(2,2) = -1.d0
      this%rotmat(10)%mat(3,3) = -1.d0
      this%rotmat(10)%name = 'C2x'
!
      this%rotmat(11)%mat(1,1) = -1.d0
      this%rotmat(11)%mat(2,2) =  1.d0
      this%rotmat(11)%mat(3,3) = -1.d0
      this%rotmat(11)%name = 'C2y'
!
      this%rotmat(12)%mat(1,1) = -1.d0
      this%rotmat(12)%mat(2,2) = -1.d0
      this%rotmat(12)%mat(3,3) =  1.d0
      this%rotmat(12)%name = 'C2z'
!
      this%rotmat(13)%mat(1,1) =  1.d0
      this%rotmat(13)%mat(2,3) =  1.d0
      this%rotmat(13)%mat(3,2) = -1.d0
      this%rotmat(13)%name = 'C4x'
!
      this%rotmat(14)%mat(1,3) = -1.d0
      this%rotmat(14)%mat(2,2) =  1.d0
      this%rotmat(14)%mat(3,1) =  1.d0
      this%rotmat(14)%name = 'C4y '
!
      this%rotmat(15)%mat(1,2) =  1.d0
      this%rotmat(15)%mat(2,1) = -1.d0
      this%rotmat(15)%mat(3,3) =  1.d0
      this%rotmat(15)%name = 'C4z'
!
      this%rotmat(16)%mat(1,1) =  1.d0
      this%rotmat(16)%mat(2,3) = -1.d0
      this%rotmat(16)%mat(3,2) =  1.d0
      this%rotmat(16)%name = 'C4x-1 '
!
      this%rotmat(17)%mat(1,3) =  1.d0
      this%rotmat(17)%mat(2,2) =  1.d0
      this%rotmat(17)%mat(3,1) = -1.d0
      this%rotmat(17)%name = 'C4y-1'
!
      this%rotmat(18)%mat(1,2) = -1.d0
      this%rotmat(18)%mat(2,1) =  1.d0
      this%rotmat(18)%mat(3,3) =  1.d0
      this%rotmat(18)%name = 'C4z-1'
!
      this%rotmat(19)%mat(1,2) =  1.d0
      this%rotmat(19)%mat(2,1) =  1.d0
      this%rotmat(19)%mat(3,3) = -1.d0
      this%rotmat(19)%name = 'C2a'
!
      this%rotmat(20)%mat(1,2) = -1.d0
      this%rotmat(20)%mat(2,1) = -1.d0
      this%rotmat(20)%mat(3,3) = -1.d0
      this%rotmat(20)%name = 'C2b'
!
      this%rotmat(21)%mat(1,3) =  1.d0
      this%rotmat(21)%mat(2,2) = -1.d0
      this%rotmat(21)%mat(3,1) =  1.d0
      this%rotmat(21)%name = 'C2c'
!
      this%rotmat(22)%mat(1,3) = -1.d0
      this%rotmat(22)%mat(2,2) = -1.d0
      this%rotmat(22)%mat(3,1) = -1.d0
      this%rotmat(22)%name = 'C2d'
!
      this%rotmat(23)%mat(1,1) = -1.d0
      this%rotmat(23)%mat(2,3) =  1.d0
      this%rotmat(23)%mat(3,2) =  1.d0
      this%rotmat(23)%name = 'C2e'
!
      this%rotmat(24)%mat(1,1) = -1.d0
      this%rotmat(24)%mat(2,3) = -1.d0
      this%rotmat(24)%mat(3,2) = -1.d0
      this%rotmat(24)%name = 'C2f'

do i1=1,24
   this%rotmat(i1+24)%mat = -this%rotmat(i1)%mat
   this%rotmat(i1+24)%name = 'I'//this%rotmat(i1)%name(1:9)
end do
!
!
!*********************************************
! Trigonal and hexagonal groups
!*********************************************
!
      this%rotmat(49)%mat(1,1) = -HALF
      this%rotmat(49)%mat(1,2) =  RTHREE
      this%rotmat(49)%mat(2,1) = -RTHREE
      this%rotmat(49)%mat(2,2) = -HALF
      this%rotmat(49)%mat(3,3) =  1.d0
      this%rotmat(49)%name = 'C3z'
!
      this%rotmat(50)%mat(1,1) = -HALF
      this%rotmat(50)%mat(1,2) = -RTHREE
      this%rotmat(50)%mat(2,1) =  RTHREE
      this%rotmat(50)%mat(2,2) = -HALF
      this%rotmat(50)%mat(3,3) =  1.d0
      this%rotmat(50)%name = 'C3z-1'
!
      this%rotmat(51)%mat(1,1) =  HALF
      this%rotmat(51)%mat(1,2) =  RTHREE
      this%rotmat(51)%mat(2,1) = -RTHREE
      this%rotmat(51)%mat(2,2) =  HALF
      this%rotmat(51)%mat(3,3) =  1.d0
      this%rotmat(51)%name = 'C6z'
!
      this%rotmat(52)%mat(1,1) =  HALF
      this%rotmat(52)%mat(1,2) = -RTHREE
      this%rotmat(52)%mat(2,1) =  RTHREE
      this%rotmat(52)%mat(2,2) =  HALF
      this%rotmat(52)%mat(3,3) =  1.d0
      this%rotmat(52)%name = 'C6z-1'
!
      this%rotmat(53)%mat(1,1) = -HALF
      this%rotmat(53)%mat(1,2) =  RTHREE
      this%rotmat(53)%mat(2,1) =  RTHREE
      this%rotmat(53)%mat(2,2) =  HALF
      this%rotmat(53)%mat(3,3) = -1.d0
      this%rotmat(53)%name = 'C2A'
!
      this%rotmat(54)%mat(1,1) = -HALF
      this%rotmat(54)%mat(1,2) = -RTHREE
      this%rotmat(54)%mat(2,1) = -RTHREE
      this%rotmat(54)%mat(2,2) =  HALF
      this%rotmat(54)%mat(3,3) = -1.d0
      this%rotmat(54)%name = 'C2B'
!
      this%rotmat(55)%mat(1,1) =  HALF
      this%rotmat(55)%mat(1,2) = -RTHREE
      this%rotmat(55)%mat(2,1) = -RTHREE
      this%rotmat(55)%mat(2,2) = -HALF
      this%rotmat(55)%mat(3,3) = -1.d0
      this%rotmat(55)%name = 'C2C'
!
      this%rotmat(56)%mat(1,1) =  HALF
      this%rotmat(56)%mat(1,2) =  RTHREE
      this%rotmat(56)%mat(2,1) =  RTHREE
      this%rotmat(56)%mat(2,2) = -HALF
      this%rotmat(56)%mat(3,3) = -1.d0
      this%rotmat(56)%name = 'C2D'

do is=1,8
   this%rotmat(56+is)%mat = -this%ROTMAT(48+is)%mat
   this%rotmat(56+is)%name = 'I'//this%ROTMAT(48+is)%name(1:9)
end do

!-----------------------------------------
end subroutine

end module
