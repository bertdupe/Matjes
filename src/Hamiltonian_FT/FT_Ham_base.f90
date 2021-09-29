module m_FT_Ham_base
use m_parameters_FT_Ham
! contains the common data to all Hamiltonian
!
!
!
implicit none

private
public FT_Ham_base

type,abstract ::    FT_Ham_base
    real(8),allocatable :: diffR(:,:)   !real-space difference of unit-cells in nm for each H_loc [3,N_H]
    integer             :: dimH=0       !dimension of Hamiltonian
    integer             :: N_H=0        !number of real-space Hamiltonians considered
    logical,private     :: set=.false.  !check if arrays have been set
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
        import FT_Ham_base,HK_inp_t,parameters_TB_IO_H
        class(FT_Ham_base),intent(inout)       :: this
        type(Hk_inp_t),intent(inout)        :: Hk_inp
        type(parameters_TB_IO_H),intent(in) :: io
    end subroutine

    subroutine int_set_work(this,work)
        import H_k_base,work_ham
        class(H_k_base),intent(inout)   :: this
        type(work_ham),intent(inout)    :: work
    end subroutine

    subroutine int_get_size_eval(this,eval_size)
        import H_k_base
        class(H_k_base),intent(inout)   :: this
        integer,intent(out)             :: eval_size
    end subroutine

    subroutine int_set_k(this,k)
        import H_k_base
        class(H_k_base),intent(inout)   :: this
        real(8),intent(in)              :: k(3)
    end subroutine

    subroutine int_eval(this,Nin,eval,Nout,work)
        import H_k_base,work_ham
        class(H_k_base),intent(inout)   :: this
        integer,intent(in)              :: Nin  !size of eigenvalue input array
        real(8),intent(inout)           :: eval(Nin)    !eigenvalue array
        integer,intent(out)             :: Nout !calculated number of eigenvalues
        type(work_ham)                      :: work !work array that should be set to correct sizes
    end subroutine

    subroutine int_evec(this,Nin,eval,evec,Nout,work)
        import H_k_base, work_ham
        class(H_k_base),intent(inout)   :: this
        integer,intent(in)                  :: Nin  !size of eigenvalue input array
        real(8),intent(inout)               :: eval(Nin)
        complex(8),intent(inout)            :: evec(this%dimH,Nin)
        integer,intent(out)                 :: Nout !calculated number of eigenvalues
        type(work_ham)                      :: work !work array that should be set to correct sizes
    end subroutine
end interface
contains

integer function  get_size_eval(this)
    class(H_k_base),intent(in)      :: this
    get_size_eval=this%dimH
end function

subroutine get_eval(this,k,Nin,eval,Nout,work)
    class(H_k_base),intent(inout)   :: this
    real(8),intent(in)              :: k(3)
    integer,intent(in)              :: Nin  !size of eigenvalue input array
    real(8),intent(inout)           :: eval(Nin)    !eigenvalue array
    integer,intent(out)             :: Nout !calculated number of eigenvalues
    type(work_ham)                  :: work !work array that should be set to correct sizes

    Call this%set_k(k)
    Call this%calc_eval(Nin,eval,Nout,work)
end subroutine

subroutine get_evec(this,k,Nin,eval,evec,Nout,work)
    class(H_k_base),intent(inout)   :: this
    real(8),intent(in)              :: k(3)
    integer,intent(in)              :: Nin  !size of eigenvalue input array
    real(8),intent(inout)           :: eval(Nin)    !eigenvalue array
    complex(8),intent(inout)        :: evec(this%dimH,Nin)    !eigenvalue array
    integer,intent(out)             :: Nout !calculated number of eigenvalues
    type(work_ham)                  :: work !work array that should be set to correct sizes

    Call this%set_k(k)
    Call this%calc_evec(Nin,eval,evec,Nout,work)
end subroutine


subroutine do_set(this,set)
    class(H_K_base),intent(inout)       :: this
    logical,intent(in)                  :: set

    this%set=set
end subroutine

logical function is_set(this)
    class(H_K_base),intent(in)  :: this
    is_set=this%set
end function

integer function get_dimH(this)
    class(H_K_base),intent(in)  :: this
    get_dimH=this%dimH
end function



end module m_FT_Ham_base
