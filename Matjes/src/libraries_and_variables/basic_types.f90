module m_basic_types

!!!!
!Operator of the simulations
!!!!

type Op_real
!     real(kind=8), pointer, contiguous :: Op_loc(:,:)
     real(kind=8), pointer :: Op_loc(:,:)
end type Op_real

type Op_Im
     complex(kind=16), pointer, contiguous :: Op_loc(:,:)
end type Op_Im

!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Basic type to store strings with different length
!
!!!!!!!!!!!!!!!!!!!!!!!!!

! varible that contains the the excitations form (sweep of EM field...)
type excitations
    real(kind=8) :: start_value(3),end_value(3)
    integer :: t_start,t_end
    character(len=30) :: name
end type excitations

! variable that contains the name of a variable
type var_name
     character(len=30) :: name
end type var_name

!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Basic type to store atoms
!
!!!!!!!!!!!!!!!!!!!!!!!!!

type atom
     character(len=30) :: name
     real(kind=8) :: moment
     real(kind=8), dimension(3) :: position
end type atom

!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Basic type to store vectors and pointers
!
!!!!!!!!!!!!!!!!!!!!!!!!!

! type symmetry operation
type symop
 real(kind=8) :: mat(3,3)
 character(len=10) :: name
end type symop

! type vector and vector pointers
type vec
     real(kind=8) :: w(3)
end type vec

type vec_dim_N
     real(kind=8), allocatable :: w(:)
end type vec_dim_N

type vec_point
     real(kind=8), pointer :: w(:)
end type vec_point

! simple pointer types
type int_pointer
     integer, pointer :: p
end type

type real_pointer
     real(kind=8), pointer :: p
end type

! this target is for the boolean variable type
type bool_var
     logical :: value
     character(len=30) :: name
end type bool_var

type vec_var
     type(real(kind=8)), dimension(3) :: value
     character(len=30) :: name
end type vec_var

type real_var
     real(kind=8) :: value
     character(len=30) :: name
end type real_var

type int_var
     integer :: value
     character(len=30) :: name
end type int_var

end module m_basic_types
