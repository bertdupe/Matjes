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
        real(8),intent(in)                 :: k(3)
    end subroutine

    subroutine int_eval(this,Nin,eval,Nout)
        import FT_Ham_base
        class(FT_Ham_base),intent(inout)      :: this
        integer,intent(in)                    :: Nin  !size of eigenvalue input array
        complex(8),intent(out)                :: eval(Nin)    !eigenvalue array
        integer,intent(out)                   :: Nout !calculated number of eigenvalues
    end subroutine

    subroutine int_evec(this,Nin,eval,evec,Nout)
        import FT_Ham_base
        class(FT_Ham_base),intent(inout)    :: this
        integer,intent(in)                  :: Nin  !size of eigenvalue input array
        complex(8),intent(out)              :: eval(Nin)
        complex(8),intent(out)              :: evec(this%dimH,Nin)
        integer,intent(out)                 :: Nout !calculated number of eigenvalues
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



end module m_FT_Ham_base
