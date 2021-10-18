module m_H_tb_base
use m_ham_init_type, only: parameters_ham_init 
private
public H_TB
type,abstract :: H_TB   !
    logical,private         :: set=.false. !has this object been set?
    !size parameters
    integer                 :: nsc=0    !size multiplicator due to BdG-trafo (set to 1 or 2)
    integer                 :: nspin=0  !size multiplicator due to spin-multiplicity (set to 1 or 2)
    integer                 :: norb=0   !size multiplicator due to number of orbitals per unit-cell
    integer                 :: ncell=0  !number of unit-cells in super-cell
    integer                 :: ndim=0   !size of basis for each unit_cell (within each BdG quadrant)
    integer                 :: dimH=0   !dimension of Hamiltonian
    !solve parameters
    integer                 :: estNe=0          !estimated number of eigenvalues  
    real(8)                 :: Ebnd(2)=0.0d0    !minimal and maximal considered eigenvalue( if such a 
    real(8)                 :: diag_acc=1d-12    ! accuracy of iterative eigenvalue solution (so far only fpm input)
contains
!    procedure(int_get_eval),deferred    :: get_eval
    procedure                           :: get_eval
    procedure(int_get_evec),deferred    :: get_evec

    procedure,NON_OVERRIDABLE       :: add
    procedure,NON_OVERRIDABLE       :: copy 
    procedure,NON_OVERRIDABLE       :: destroy
    procedure,NON_OVERRIDABLE       :: init_base
    procedure,NON_OVERRIDABLE       :: init_otherH
    procedure,NON_OVERRIDABLE       :: get_hinit

    procedure(int_mult),deferred    :: mult_r

    procedure(int_init_connect),deferred    :: init_connect
    procedure(int_init_coo),deferred        :: init_coo
    procedure(int_mv),deferred              :: mv

    procedure(int_destroy),deferred :: destroy_child
    procedure(int_add_H),deferred   :: add_child
    procedure(int_copy),deferred    :: copy_child
    !misc.
    procedure,NON_OVERRIDABLE  :: is_set
    procedure,NON_OVERRIDABLE  :: set_prepared
end type

abstract interface

    subroutine int_mult(this,vec,res,alpha,beta)
        import H_tb
        class(H_tb),intent(in)          :: this
        complex(8),intent(in   )        :: vec(this%dimH)
        complex(8),intent(inout)        :: res(this%dimH)
        complex(8),intent(in),optional  :: alpha
        complex(8),intent(in),optional  :: beta
    end subroutine

    subroutine int_mv(this,Hout)
        import H_TB
        class(H_tb),intent(inout)  ::  this
        class(H_tb),intent(inout)  ::  Hout
    end subroutine


    subroutine int_get_evec(this,eval,evec)
        import H_TB
        class(H_tb),intent(in)  ::  this
        real(8),allocatable,intent(out)     ::  eval(:)
        complex(8),allocatable,intent(out)  ::  evec(:,:)
    end subroutine


    subroutine int_init_connect(this,connect,Hval,Hval_ind,io,ind_offset)
        import H_TB, parameters_ham_init
        class(H_tb),intent(inout)           :: this
        complex(8),intent(in)               :: Hval(:) 
        integer,intent(in)                  :: Hval_ind(:,:)
        integer,intent(in)                  :: connect(:,:) 
        type(parameters_ham_init),intent(in):: io
        integer,intent(in),optional         :: ind_offset(2)
    end subroutine

    subroutine int_init_coo(this,val,row,col,io,ind_offset)
        import H_TB, parameters_ham_init
        class(H_tb),intent(inout)           :: this
        complex(8),intent(in)               :: val(:) 
        integer,intent(in)                  :: row(size(val)),col(size(val))
        type(parameters_ham_init),intent(in):: io
        integer,intent(in),optional         :: ind_offset(2)
    end subroutine

    subroutine int_destroy(this)
        import H_TB
        class(H_TB),intent(inout)  :: this
    end subroutine

    subroutine int_add_H(this,H_in)
        import H_TB
        class(H_TB),intent(inout)    :: this
        class(H_TB),intent(in)       :: H_in
    end subroutine

    subroutine int_copy(this,Hout)
        import H_TB
        class(H_TB),intent(in)         :: this
        class(H_TB),intent(inout)      :: Hout
    end subroutine
end interface

contains

subroutine init_otherH(this,H_in)
    class(H_TB),intent(inout)    :: this
    class(H_TB),intent(in)       :: h_in

    this%dimH    =H_in%dimH 
    this%nsc     =H_in%nsc  
    this%nspin   =H_in%nspin
    this%norb    =H_in%norb 
    this%ncell   =H_in%ncell
    this%ndim    =H_in%ndim 
    this%Ebnd    =H_in%Ebnd 
    this%estNe   =H_in%estNe 
    this%diag_acc=H_in%diag_acc 
    this%set=.true.
end subroutine

subroutine init_base(this,io)
    class(H_TB),intent(inout)           :: this
    type(parameters_ham_init),intent(in) :: io

    this%nsc     =io%nsc  
    this%nspin   =io%nspin  
    this%norb    =io%norb  
    this%ncell   =io%ncell  

    this%Ebnd    =io%Ebnd
    this%estNe   =io%estNe
    this%diag_acc=io%diag_acc

    this%ndim    =this%norb*this%nspin  
    this%dimH    =this%norb*this%nspin*this%nsc*this%ncell
    this%set =.true.
end subroutine

subroutine get_hinit(this,hinit)
    class(H_TB),intent(in)                  :: this
    type(parameters_ham_init),intent(inout) :: hinit

    hinit%nsc     =this%nsc  
    hinit%nspin   =this%nspin  
    hinit%norb    =this%norb  
    hinit%ncell   =this%ncell  

    hinit%Ebnd    =this%Ebnd
    hinit%estNe   =this%estNe
    hinit%diag_acc=this%diag_acc
end subroutine



function is_set(this) result(l)
    class(H_TB),intent(in)      ::  this
    logical                     ::  l

    l=this%set
end function

subroutine set_prepared(this,l)
    class(H_TB),intent(inout)   ::  this
    logical,intent(in)          ::  l

    this%set=l
end subroutine


recursive subroutine copy(this,Hout)
    !copy this to Hout
    class(H_TB),intent(in)         :: this
    class(H_TB),intent(inout)      :: Hout

    if(this%is_set())then
        Call Hout%destroy()

        Hout%dimH    =this%dimH
        Hout%nsc     =this%nsc 
        Hout%nspin   =this%nspin 
        Hout%norb    =this%norb 
        Hout%ncell   =this%ncell 
        Hout%ndim    =this%ndim 
        Hout%Ebnd    =this%Ebnd 
        Hout%estNe   =this%estNe 
        Hout%diag_acc=this%diag_acc 

        Call this%copy_child(Hout)
        Call Hout%set_prepared(.true.)
    else
        ERROR STOP "cannot copy H since source is not set"
    endif
end subroutine


recursive subroutine add(this,H_in)
    class(H_TB),intent(inout)    :: this
    class(H_TB),intent(in)       :: H_in

    if(this%is_set())then
        if(.not.H_in%is_set())     ERROR STOP "CANNOT ADD hamiltonians as H_in is not set"
        if(this%dimH/=H_in%dimH)   ERROR STOP "CANNOT ADD Hamiltonians with different Hamiltonian dimensions"
        if(this%nspin/=H_in%nspin) ERROR STOP "CANNOT ADD Hamiltonians with different number of spins"
        if(this%norb/=H_in%norb)   ERROR STOP "CANNOT ADD Hamiltonians with different number of orbitals"
        if(this%ncell/=H_in%ncell) ERROR STOP "CANNOT ADD Hamiltonians with different number of unit-cells"
        Call this%add_child(H_in)
    else
        Call H_in%copy(this)
    endif
end subroutine

subroutine destroy(this)
    !Destroys Hamiltonian
    class(H_TB),intent(inout)    :: this

    this%dimH=0
    this%nsc=0 
    this%nspin=0 
    this%norb=0 
    this%ncell=0 
    this%ndim=0 
    this%estNe=0
    this%diag_acc=0
    this%Ebnd=0.0d0 
    Call this%destroy_child()
    this%set=.false.
end subroutine

subroutine get_eval(this,eval)
    use, intrinsic :: iso_fortran_env, only : error_unit
    !get eigenvalues from eigenvectors, overwrite if diagonalization methods allows for calculating eigenvalues only
    class(H_TB),intent(in)          ::  this
    real(8),intent(out),allocatable ::  eval(:)
    complex(8),allocatable          ::  evec(:,:)
    logical,save                    ::  said=.false.

    if(.not.said)then
        write(error_unit,'(2/A)') "WARNING, requested eigenvalues only but chosen diagonalization routine only supports calculation including the eigenvectors"
        write(error_unit,'(A/)')  "Another diagonalization routine might be more efficient"
        said=.true.
    endif
    Call this%get_evec(eval,evec)
    deallocate(evec)
end subroutine



end module
