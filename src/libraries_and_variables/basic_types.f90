module m_basic_types

private :: atom_bcast
!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Basic type to store strings with different length
!
!!!!!!!!!!!!!!!!!!!!!!!!!

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
     real(8), dimension(3) :: position
     real(8) :: moment=0.0d0
     real(8) :: charge=0.0d0
     contains
     procedure :: bcast => atom_bcast
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

!! type vector and vector pointers
!type vec
!     real(kind=8) :: w(3)
!end type vec

type vec_dim_N
     real(kind=8), allocatable :: w(:)
end type vec_dim_N

!type vec_point
!     real(kind=8), pointer :: w(:)
!end type vec_point

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
     contains
     procedure :: bcast => bool_var_bcast
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
contains 

subroutine atom_bcast(this,comm)
    use mpi_basic                
    class(atom),intent(inout)    ::  this
    type(mpi_type),intent(in)       ::  comm

#ifdef CPP_MPI
    integer     :: ierr
    Call MPI_Bcast(this%name, 30, MPI_CHARACTER, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%position, 3, MPI_REAL8, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%moment, 1, MPI_REAL8, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%charge, 1, MPI_REAL8, comm%mas, comm%com,ierr)
    !could be done more elegantly with custom MPI_type
#else
    continue
#endif
end subroutine


subroutine bool_var_bcast(this,comm)
    use mpi_basic                
    class(bool_var),intent(inout)   ::  this
    type(mpi_type),intent(in)       ::  comm

#ifdef CPP_MPI
    integer     :: ierr
    Call MPI_Bcast(this%name, 30, MPI_CHARACTER, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%value, 1, MPI_LOGICAL, comm%mas, comm%com,ierr)
    !could be done more elegantly with custom MPI_type
#else
    continue
#endif
end subroutine

end module m_basic_types
