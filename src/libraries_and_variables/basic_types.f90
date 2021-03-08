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

type atomtype
    !type for reading the input
    character(len=30) :: name   !some reference of the name of the Atom type(can be used to specify atomtype in unitcell)
    real(8) :: moment=0.0d0 !magnetic moment
    real(8) :: charge=0.0d0 !effective charge 
    real(8) :: mass=0.0d0   !effective mass for phonons
    logical :: use_ph=.false.   !use for phonon calculation
    integer :: orbitals=0   !number of tight-binding orbitals per spin
    contains
    procedure :: atomtype_write
    generic :: write(formatted) => atomtype_write   !this is rather nonsensical and should rather be used to write it the atom-type in some internal format
end type

type atom
    !type which contains all essential information about the an atom in the unit-cell
    real(8), dimension(3)   :: position
    integer ::  type_id=0    !index of the atom type(should be 1..N_atomtypes) (can be used to specify atomtype in unitcell, and check if same atom type is used)
    character(len=30) :: name   !some reference of the name of the Atom type(can be used to specify atomtype in unitcell)
    real(8) :: moment=0.0d0 !magnetic moment
    real(8) :: charge=0.0d0 !effective charge 
    real(8) :: mass=0.0d0   !effective mass for phonons
    logical :: use_ph=.false.   !use for phonon calculation
    integer :: orbitals=0   !number of tight-binding orbitals per spin
    contains
    procedure :: bcast => atom_bcast
    procedure :: set_attype => atom_set_attype
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
    Call MPI_Bcast(this%position, 3             , MPI_DOUBLE_PRECISION, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%type_id , 1             , MPI_INTEGER         , comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%name    , len(this%name), MPI_CHARACTER       , comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%moment  , 1             , MPI_DOUBLE_PRECISION, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%charge  , 1             , MPI_DOUBLE_PRECISION, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%mass    , 1             , MPI_DOUBLE_PRECISION, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%use_ph  , 1             , MPI_LOGICAL         , comm%mas, comm%com,ierr)
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

subroutine atomtype_write(attype, unit, iotype, v_list, iostat, iomsg)
  class(atomtype), intent(in)   :: attype
  integer, intent(in)           :: unit
  character(*), intent(in)      :: iotype
  integer, intent(in)           :: v_list(:)
  integer, intent(out)          :: iostat
  character(*), intent(inout)   :: iomsg

  write (unit, "(2A,/3(A,F8.4/),A,L3/,A,I3)", IOSTAT=iostat, IOMSG=iomsg)  &
      "atom type:  name     : ",trim(attype%name),&
      "            moment   : ",attype%moment,&
      "            charge   : ",attype%charge,&
      "            mass     : ",attype%mass,&
      "            use ph   : ",attype%use_ph,&
      "            #orbitals: ",attype%orbitals
end subroutine

subroutine atom_set_attype(this,attype)
    class(atom),intent(inout)   ::  this
    class(atomtype),intent(in)  ::  attype

    this%name=attype%name
    this%moment=attype%moment
    this%charge=attype%charge
    this%mass=attype%mass
    this%use_ph=attype%use_ph
    this%orbitals=attype%orbitals
end subroutine
end module m_basic_types
