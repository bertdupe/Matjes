module m_FT_Ham_base
use m_parameters_FT_Ham
use m_FT_Ham_coo_rtok_base
! contains the common data to all Hamiltonian
!
!
!
implicit none

private
public FT_Ham_base

type,abstract ::    FT_Ham_base
    integer                      :: dimH=0       !dimension of Hamiltonian
    integer                      :: N_H=0        !number of real-space Hamiltonians considered
    logical,private              :: set=.false.  !check if arrays have been set
    type(parameters_FT_HAM_IO)   :: io_H         ! parameters for the diagonalization
contains
    procedure(int_init),deferred        :: init        !initializes this type
    procedure(int_set_k),deferred       :: set_k       !sets the internal Hamiltonian to a k-point
    procedure(int_set_work),deferred    :: set_work    !sets the work arrays to the correct size
    procedure(int_eval),deferred        :: calc_eval   !subroutine to calculate the eigenvalues using the internally set matrix
    procedure(int_evec),deferred        :: calc_evec   !subroutine to calculate the eigenvectors and eigenvalues using the internally set matrix

    procedure                   :: get_size_eval
    procedure,NON_OVERRIDABLE   :: get_dimH
    procedure,NON_OVERRIDABLE   :: get_eval     !get eigenvalues  for a given k
    procedure,NON_OVERRIDABLE   :: get_evec     !get eigenvectors for a given k
    procedure,NON_OVERRIDABLE   :: do_set       !set set logical
    procedure,NON_OVERRIDABLE   :: is_set       !returns set logical

    !MPI-stuff
    procedure,NON_OVERRIDABLE               :: send_base
    procedure,NON_OVERRIDABLE               :: recv_base
    procedure,NON_OVERRIDABLE               :: bcast_base
    procedure(int_send),deferred            :: send
    procedure(int_recv),deferred            :: recv
    procedure(int_mpicom),deferred          :: distribute
    procedure(int_mpicom),deferred          :: bcast

end type


abstract interface
    subroutine int_init(this,Hk_inp,io)
        import FT_Ham_base,parameters_FT_HAM_IO,H_inp_real_to_k
        class(FT_Ham_base),intent(inout)         :: this
        type(H_inp_real_to_k),intent(in)         :: Hk_inp(:)
        type(parameters_FT_HAM_IO),intent(in)    :: io
    end subroutine

    subroutine int_set_work(this,eval,evec)
        import FT_Ham_base
        class(FT_Ham_base),intent(inout)   :: this
        complex(8),allocatable,intent(out) :: eval(:)
        complex(8),allocatable,intent(out) :: evec(:,:)
    end subroutine

    subroutine int_get_size_eval(this,eval_size)
        import FT_Ham_base
        class(FT_Ham_base),intent(inout)   :: this
        integer,intent(out)                :: eval_size
    end subroutine

    subroutine int_set_k(this,Hk_inp,k)
        import FT_Ham_base,H_inp_real_to_k
        class(FT_Ham_base),intent(inout)   :: this
        type(H_inp_real_to_k),intent(in)   :: Hk_inp(:)
        real(8),intent(in)                 :: k(:)
    end subroutine

    subroutine int_eval(this,Nin,eval,Nout)
        import FT_Ham_base
        class(FT_Ham_base),intent(inout)      :: this
        integer,intent(in)                    :: Nin  !size of eigenvalue input array
        complex(8),intent(out)                :: eval(:)    !eigenvalue array
        integer,intent(out)                   :: Nout !calculated number of eigenvalues
    end subroutine

    subroutine int_evec(this,Nin,eval,evec,Nout)
        import FT_Ham_base
        class(FT_Ham_base),intent(inout)    :: this
        integer,intent(in)                  :: Nin  !size of eigenvalue input array
        complex(8),intent(out)              :: eval(:)
        complex(8),intent(out)              :: evec(:,:)
        integer,intent(out)                 :: Nout !calculated number of eigenvalues
    end subroutine

    subroutine int_mpicom(this,comm)
        use mpi_basic
        import FT_Ham_base
        class(FT_Ham_base),intent(inout)   ::  this
        type(mpi_type),intent(in)          ::  comm
    end subroutine

    subroutine int_send(this,ithread,tag,com)
        import FT_Ham_base
        class(FT_Ham_base),intent(in)      :: this
        integer,intent(in)                 :: ithread
        integer,intent(in)                 :: tag
        integer,intent(in)                 :: com
    end subroutine

    subroutine int_recv(this,ithread,tag,com)
        import FT_Ham_base
        class(FT_Ham_base),intent(inout)   :: this
        integer,intent(in)                 :: ithread
        integer,intent(in)                 :: tag
        integer,intent(in)                 :: com
    end subroutine
end interface
contains

integer function  get_size_eval(this)
    class(FT_Ham_base),intent(in)      :: this
    get_size_eval=this%dimH
end function

subroutine get_eval(this,Nin,eval,Nout)
    class(FT_Ham_base),intent(inout)   :: this
    integer,intent(in)                 :: Nin  !size of eigenvalue input array
    complex(8),intent(inout)           :: eval(Nin)    !eigenvalue array
    integer,intent(out)                :: Nout !calculated number of eigenvalues

    Call this%calc_eval(Nin,eval,Nout)
end subroutine

subroutine get_evec(this,Nin,eval,evec,Nout)
    class(FT_Ham_base),intent(inout)   :: this
    integer,intent(in)                 :: Nin  !size of eigenvalue input array
    complex(8),intent(inout)           :: eval(Nin)    !eigenvalue array
    complex(8),intent(inout)           :: evec(this%dimH,Nin)    !eigenvalue array
    integer,intent(out)                :: Nout !calculated number of eigenvalues

    Call this%calc_evec(Nin,eval,evec,Nout)
end subroutine


subroutine do_set(this,set)
    class(FT_Ham_base),intent(inout)       :: this
    logical,intent(in)                     :: set

    this%set=set
end subroutine

logical function is_set(this)
    class(FT_Ham_base),intent(in)  :: this
    is_set=this%set
end function

integer function get_dimH(this)
    class(FT_Ham_base),intent(in)  :: this
    get_dimH=this%dimH
end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!            MPI ROUTINES           !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine bcast_base(this,comm)
        use mpi_basic
        class(FT_Ham_base),intent(inout)        ::  this
        type(mpi_type),intent(in)               ::  comm
#ifdef CPP_MPI
        integer     :: ierr
        integer     :: N(2)
        integer     :: mode_ident(2)

        Call MPI_Bcast(this%dimH         ,   1, MPI_INTEGER  , comm%mas, comm%com,ierr)
        Call MPI_Bcast(this%N_H          ,   1, MPI_INTEGER  , comm%mas, comm%com,ierr)
        Call MPI_Bcast(this%set          ,   1, MPI_LOGICAL  , comm%mas, comm%com,ierr)

#else
        continue
#endif
    end subroutine

subroutine send_base(this,ithread,tag,com)
    use mpi_basic
    class(FT_Ham_base),intent(in)   :: this
    integer,intent(in)              :: ithread
    integer,intent(in)              :: tag
    integer,intent(in)              :: com

#ifdef CPP_MPI
    integer     :: ierr

    if(.not.this%is_set()) ERROR STOP "CANNOT SEND HAMILTONIAN WHEN THE ORIGIN IS NOT SET"


    Call MPI_SEND(this%dimH,          1,               MPI_INTEGER,   ithread, tag, com, ierr)
    Call MPI_SEND(this%N_H,           1,               MPI_INTEGER,   ithread, tag, com, ierr)
    Call MPI_SEND(this%set,           1,               MPI_LOGICAL,   ithread, tag, com, ierr)

#else
    continue
#endif
end subroutine

subroutine recv_base(this,ithread,tag,com)
    use mpi_basic
    class(FT_Ham_base),intent(inout)  :: this
    integer,intent(in)                :: ithread
    integer,intent(in)                :: tag
    integer,intent(in)                :: com

#ifdef CPP_MPI
    integer     :: ierr
    integer     :: stat(MPI_STATUS_SIZE)

    if(this%is_set()) ERROR STOP "CANNOT RECV HAMILTONIAN WHEN THE ORIGIN IS ALREADY SET"

    Call MPI_RECV(this%dimH,          1,               MPI_INTEGER,   ithread, tag, com, stat, ierr)
    Call MPI_RECV(this%N_H,           1,               MPI_INTEGER,   ithread, tag, com, stat, ierr)
    Call MPI_RECV(this%set,           1,               MPI_LOGICAL,   ithread, tag, com, stat, ierr)

#else
    continue
#endif
end subroutine

end module m_FT_Ham_base
