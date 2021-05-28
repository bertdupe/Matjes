module m_H_type
!module containing the basic polymorphic Hamiltonian class t_H_base
use m_derived_types, only : lattice,number_different_order_parameters
use m_mode_public, only: F_mode, get_mode_ident, set_mode_ident
use m_work_ham_single, only:  work_ham_single, work_mode
use,intrinsic :: ISO_FORTRAN_ENV, only: error_unit
implicit none

private
public :: t_H_base
public :: len_desc

integer,parameter       ::  len_desc=100    !length of the character describing the Hamiltonian

type,abstract :: t_H_base
    !abstrace base class for the Hamiltonian 
    integer                     :: dimH(2)=0        !dimension of Hamiltonian
    integer,allocatable         :: op_l(:),op_r(:)  !operator indices which have to be combined to get the space of the left/right side of the Hamiltonian
    integer                     :: dim_mode(2)=0    !size of left/right mode after multiplying out order parameters (size per lattice point)
    class(F_mode),allocatable   :: mode_l,mode_r    !class to get the left/right mode on which the Hamiltonian acts
    logical,private             :: set=.false.      !has this object been set?
    character(len=len_desc)     :: desc=""          !description of the Hamiltonian term, only used for user information and should be set manually 
    integer                     :: mult_M_single=0  !factor necessary to calculate energy change correctly when only evaluating single sites

    !variables necessary to sparse matrix ind evaluation    (is not set correctly for all Hamiltonians)
    integer                     :: col_max=0 !maximal number of entries for left  multiplication
    integer                     :: row_max=0 !maximal number of entries for right multiplication
    !sizes of single evaluation for respective order parameter
    integer                     :: dim_l_single(number_different_order_parameters)=0 !dimension for inner work size array of set order on left side(set_work_size_single)
    integer                     :: dim_r_single(number_different_order_parameters)=0 !dimension for inner work size array of set order (set_work_size_single)
contains

    !Hamiltonian initialization routines
    procedure(int_init_H_connect),deferred          :: init_connect         !based on Hamiltonian in coo format input (arrays), only rank_2, connection input
    procedure(int_init_H_mult_connect_2),deferred   :: init_mult_connect_2  !based on Hamiltonian in coo format input (arrays), for rank >2, but only at 2 sites at once, connection input
    procedure(int_init_H_coo),deferred              :: init_coo             !based on Hamiltonian in coo format input (arrays), the coo data-arrays are directly moved into the type
    procedure(int_destroy),deferred                 :: optimize             !calls possible optimization routines rearanging internal Hamiltonian layout

    !routines acting on the entire Hamiltonian
    procedure ,NON_OVERRIDABLE              :: eval_all                 !evaluates energy of full lattice
    procedure(int_mult),deferred            :: mult_r,mult_l            !multipy out with left/right side
    procedure ,NON_OVERRIDABLE              :: energy_dist              !get energy distribution per site

    !routines acting on parts of the Hamiltonian
    procedure                               :: mult_l_disc,   mult_r_disc

    !routines setting work array sizes
    procedure                               :: set_work_single      !sets the necessary sizes for the work arrays
    procedure                               :: set_work_mode        !sets the necessary sizes for higher rank modes

    !Utility functions
    procedure,NON_OVERRIDABLE               :: finish_setup_base
    procedure,NON_OVERRIDABLE               :: destroy
    procedure,NON_OVERRIDABLE               :: copy
    procedure,NON_OVERRIDABLE               :: add
    procedure,NON_OVERRIDABLE               :: init_base
    procedure,NON_OVERRIDABLE               :: init_otherH
    procedure(int_finish_setup),deferred    :: finish_setup
    procedure(int_add_H),deferred           :: add_child
    procedure(int_destroy),deferred         :: destroy_child
    procedure(int_copy),deferred            :: copy_child

    !MPI-stuff
    procedure,NON_OVERRIDABLE               :: send_base
    procedure,NON_OVERRIDABLE               :: recv_base
    procedure,NON_OVERRIDABLE               :: bcast_base
    procedure(int_send),deferred            :: send
    procedure(int_recv),deferred            :: recv
    procedure(int_mpicom),deferred          :: distribute
    procedure(int_mpicom),deferred          :: bcast

    !misc.
    procedure,NON_OVERRIDABLE  :: is_set
    procedure,NON_OVERRIDABLE  :: set_prepared
    procedure,NON_OVERRIDABLE  :: same_space

    !finalize routine? (might be risky with Hamiltonian references that have been passed around)
end type


abstract interface
    subroutine int_finish_setup(this)
        import t_H_base
        class(t_H_base),intent(inout)    :: this
    end subroutine

    subroutine int_mult(this,lat,res,work,alpha,beta)
        import t_H_base,lattice, work_mode
        class(t_H_base),intent(in)      :: this
        type(lattice),intent(in)        :: lat
        real(8),intent(inout)           :: res(:)
        type(work_mode),intent(inout)   :: work
        real(8),intent(in),optional     :: alpha
        real(8),intent(in),optional     :: beta
    end subroutine

    subroutine int_destroy(this)
        import t_H_base
        class(t_H_base),intent(inout)  :: this
    end subroutine

    subroutine int_init_H_mult_connect_2(this,connect,Hval,Hval_ind,op_l,op_r,lat,mult_M_single,dim_mode_in)
        import t_H_base,lattice
        class(t_H_base),intent(inout)   :: this
        type(lattice),intent(in)        :: lat
        character(*),intent(in)         :: op_l
        character(*),intent(in)         :: op_r
        real(8),intent(in)              :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
        integer,intent(in)              :: Hval_ind(:,:)
        integer,intent(in)              :: connect(:,:)
        integer,intent(in)              :: mult_M_single
        integer,intent(in),optional     :: dim_mode_in(2)   !optional way of putting in dim_mode directly (mainly for custom(not fully unfolded)rankN tensors)
    end subroutine

    subroutine int_init_H_coo(this,rowind,colind,val,dim_mode,op_l,op_r,lat,mult_M_single)
        import t_H_base,lattice
        class(t_H_base),intent(inout)       :: this
        real(8),allocatable,intent(inout)   :: val(:)
        integer,allocatable,intent(inout)   :: rowind(:),colind(:)
        integer,intent(in)                  :: dim_mode(2)
        character(len=*),intent(in)         :: op_l         !which order parameters are used at left  side of local Hamiltonian-matrix
        character(len=*),intent(in)         :: op_r         !which order parameters are used at right side of local Hamiltonian-matrix
        type(lattice),intent(in)            :: lat
        integer,intent(in)                  :: mult_M_single !gives the multiple with which the energy_single calculation has to be multiplied (1 for on-site terms, 2 for eg. magnetic exchange)
    end subroutine

    subroutine int_init_H_connect(this,connect,Hval,Hval_ind,order,lat,mult_M_single)
        import t_H_base,lattice
        class(t_H_base),intent(inout)   :: this
        type(lattice),intent(in)        :: lat
        character(2),intent(in)         :: order
        real(8),intent(in)              :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
        integer,intent(in)              :: Hval_ind(:,:)
        integer,intent(in)              :: connect(:,:)
        integer,intent(in)              :: mult_M_single
    end subroutine

    subroutine int_add_H(this,H_in)
        import t_H_base
        class(t_H_base),intent(inout)   :: this
        class(t_H_base),intent(in)      :: H_in
    end subroutine

    subroutine int_copy(this,Hout)
        import t_H_base
        class(t_H_base),intent(in)      :: this
        class(t_H_base),intent(inout)   :: Hout
    end subroutine

    subroutine int_mpicom(this,comm)
        use mpi_basic                
        import t_H_base
        class(t_H_base),intent(inout)   ::  this
        type(mpi_type),intent(in)       ::  comm
    end subroutine

    subroutine int_send(this,ithread,tag,com)
        import t_H_base
        class(t_H_base),intent(in)      :: this
        integer,intent(in)              :: ithread
        integer,intent(in)              :: tag
        integer,intent(in)              :: com
    end subroutine

    subroutine int_recv(this,ithread,tag,com)
        import t_H_base
        class(t_H_base),intent(inout)   :: this
        integer,intent(in)              :: ithread
        integer,intent(in)              :: tag
        integer,intent(in)              :: com
    end subroutine
end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!                    MOST USEFULL ROUTINES                        !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine eval_all(this,E,lat,work)
        USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_DOUBLE
        use, intrinsic :: iso_fortran_env, only : error_unit
        !evaluates the energie summed for all sites
        class(t_H_base),intent(in)      :: this
        type(lattice), intent(in)       :: lat
        type(work_mode),intent(inout)   :: work
        real(8), intent(out)            :: E
        ! internal
        real(C_DOUBLE)              :: tmp(this%dimH(1))
        real(8),pointer             :: modes_l(:)
    
        if(.not.allocated(this%mode_l).or..not.allocated(this%mode_r))then
            write(error_unit,'(2/2A)') "Failed to evaluate Hamiltonian: ", this%desc
            ERROR STOP "IMPLEMENT mode_l/mode_r for all Hamiltonians"
        endif

        !try to also save tmp in work (difficult because of overlap
        Call this%mult_r(lat,tmp,work)
        Call this%mode_l%get_mode(lat,modes_l,work)
        E=dot_product(modes_l,tmp)
    end subroutine 

    subroutine energy_dist(this,lat,order,work,E)
        use m_type_lattice, only: dim_modes_inner
        !get the energy values per site for a given order
        !might be wrong for non-symmetrical Hamiltonians
        class(t_H_base),intent(in)       :: this
        type(lattice), intent(in)       :: lat
        integer,intent(in)              :: order
        type(work_mode),intent(inout)   :: work
        real(8),intent(out)             :: E(lat%Ncell*lat%site_per_cell(order))
        ! internal
        real(8),pointer             :: modes(:)
        real(8),allocatable         :: tmp(:)
        real(8),allocatable         :: red(:)

        if(any(this%op_l==order))then
            allocate(tmp(this%dimH(1)))
            Call this%mult_r(lat,tmp,work)
            Call this%mode_l%get_mode(lat,modes,work)
            tmp=modes*tmp
            nullify(modes)
            if(size(this%op_l)>1)then
                allocate(red(lat%Ncell*lat%dim_modes(order)))
                Call this%mode_l%mode_reduce(lat,tmp,order,red)
                Call get_sum(dim_modes_inner(order),lat%Ncell*lat%site_per_cell(order),red,E)
                deallocate(red)
            else
                Call get_sum(dim_modes_inner(order),lat%Ncell*lat%site_per_cell(order),tmp,E)
            endif
            deallocate(tmp)
        elseif(any(this%op_r==order))then
            allocate(tmp(this%dimH(2)))
            Call this%mult_l(lat,tmp,work)
            Call this%mode_r%get_mode(lat,modes,work)
            tmp=modes*tmp
            nullify(modes)
            if(size(this%op_r)>1)then
                allocate(red(lat%Ncell*lat%dim_modes(order)))
                Call this%mode_r%mode_reduce(lat,tmp,order,red)
                Call get_sum(dim_modes_inner(order),lat%Ncell*lat%site_per_cell(order),red,E)
                deallocate(red)
            else
                Call get_sum(dim_modes_inner(order),lat%Ncell*lat%site_per_cell(order),tmp,E)
            endif
            deallocate(tmp)
        else
            E=0.0d0
        endif
    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!                    BASE ROUTINES                                !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine bcast_base(this,comm)
        use mpi_basic                
        class(t_H_base),intent(inout)        ::  this
        type(mpi_type),intent(in)       ::  comm
#ifdef CPP_MPI
        integer     :: ierr
        integer     :: N(2)
        integer     :: mode_ident(2)

        Call MPI_Bcast(this%dimH         ,   2, MPI_INTEGER  , comm%mas, comm%com,ierr)
        Call MPI_Bcast(this%dim_mode     ,   2, MPI_INTEGER  , comm%mas, comm%com,ierr)
        Call MPI_Bcast(this%mult_M_single,   2, MPI_INTEGER  , comm%mas, comm%com,ierr)
        Call MPI_Bcast(this%set          ,   1, MPI_LOGICAL  , comm%mas, comm%com,ierr)
        Call MPI_Bcast(this%desc         , 100, MPI_CHARACTER, comm%mas, comm%com,ierr)

        if(comm%ismas)then
            N(1)=size(this%op_l)
            N(2)=size(this%op_r)
        endif
        Call MPI_Bcast(N, 2, MPI_INTEGER, comm%mas, comm%com,ierr)
        if(.not.comm%ismas)then
            allocate(this%op_l(N(1)),this%op_r(N(2)))
        endif
        Call MPI_Bcast(this%op_l,N(1), MPI_INTEGER, comm%mas, comm%com,ierr)
        Call MPI_Bcast(this%op_r,N(2), MPI_INTEGER, comm%mas, comm%com,ierr)

        !bcast modes
        if(comm%ismas)then
            Call get_mode_ident(this%mode_l,mode_ident(1))
            Call get_mode_ident(this%mode_r,mode_ident(2))
        endif
        Call MPI_Bcast(mode_ident, 2, MPI_INTEGER, comm%mas, comm%com,ierr)
        if(.not.comm%ismas)then
            Call set_mode_ident(this%mode_l,mode_ident(1))
            Call set_mode_ident(this%mode_r,mode_ident(2))
        endif
        Call this%mode_l%bcast(comm)
        Call this%mode_r%bcast(comm)
        Call MPI_Bcast(this%dim_l_mode, number_different_order_parameters, MPI_INTEGER  , comm%mas, comm%com,ierr)
        Call MPI_Bcast(this%dim_r_mode, number_different_order_parameters, MPI_INTEGER  , comm%mas, comm%com,ierr)
        Call this%finish_setup()
#else
        continue
#endif
    end subroutine

    subroutine copy(this,Hout)
        !copy this to Hout
        class(t_H_base),intent(in)         :: this
        class(t_H_base),intent(inout)      :: Hout

        if(this%is_set())then
            if(.not.allocated(this%op_l).or..not.allocated(this%op_r))&
                & STOP "Cannot copy H since base op_l/op_r are not allocated"
            Call Hout%destroy()
            allocate(Hout%op_l,source=this%op_l)
            allocate(Hout%op_r,source=this%op_r)
            Hout%dimH=this%dimH
            Hout%dim_mode=this%dim_mode
            Call this%copy_child(Hout)
            Call Hout%set_prepared(.true.)
            Hout%desc=this%desc
            Hout%mult_M_single=this%mult_M_single
            if(allocated(this%mode_l)) Call this%mode_l%copy(Hout%mode_l)
            if(allocated(this%mode_r)) Call this%mode_r%copy(Hout%mode_r)
            Hout%dim_l_single=this%dim_l_single
            Hout%dim_r_single=this%dim_r_single
            if(allocated(this%mode_l).and.allocated(this%mode_r)) Call Hout%finish_setup()
        else
            STOP "cannot copy H since source is not set"
        endif
    end subroutine

    subroutine add(this,H_in)
        use m_type_lattice, only : op_int_to_abbrev
        !add H_in to this Hamiltonian, or set this to H_in is this is not set
        class(t_H_base),intent(inout)    :: this
        class(t_H_base),intent(in)       :: H_in

        if(this%is_set())then
            !check that it is possible to add Hamiltonians
            if(.not.allocated(H_in%op_l).or..not.allocated(H_in%op_r)) &
                & STOP "cannot add H and its op_l/op_r-arguments are not allocated"
            if(.not.all(this%op_l==H_in%op_l)) ERROR STOP "CANNOT ADD hamiltonians with different op_l"
            if(.not.all(this%op_r==H_in%op_r)) ERROR STOP "CANNOT ADD hamiltonians with different op_r"
            if(any(this%dimH/=H_in%dimH)) ERROR STOP "Trying to add Hamiltonians with different Hamiltonian dimensions"
            Call check_mode_same(this%mode_l,H_in%mode_l)
            Call check_mode_same(this%mode_r,H_in%mode_r)

            Call this%add_child(H_in)
            this%desc="sum in "//trim(op_int_to_abbrev(this%op_l))//trim(op_int_to_abbrev(this%op_r))//"-space"
            if(this%mult_M_single/=H_in%mult_M_single) this%mult_M_single=0 !set single multiply prefactor to 0 in order to notice if this is evaluated incorrectly
        else
            Call H_in%copy(this)
        endif
    end subroutine

    subroutine destroy(this)
        !Destroys Hamiltonian
        class(t_H_base),intent(inout)    :: this

        Call this%destroy_child()
        this%dimH=0
        if(allocated(this%op_l)) deallocate(this%op_l)
        if(allocated(this%op_r)) deallocate(this%op_r)
        this%dim_mode=0
        Call this%set_prepared(.false.)
        this%desc=""
        this%mult_M_single=0
        if(allocated(this%mode_l))then
            Call this%mode_l%destroy()
            deallocate(this%mode_l)
        endif
        if(allocated(this%mode_r))then
            Call this%mode_r%destroy()
            deallocate(this%mode_r)
        endif
        this%dim_l_single=0
        this%dim_r_single=0
    end subroutine

    subroutine init_base(this,lat,op_l,op_r)
        class(t_H_base),intent(inout)    :: this
        class(lattice),intent(in)   :: lat
        !might make sense to change op_l, op_r to characters
        integer,intent(in)          :: op_l(:),op_r(:)
        integer     ::  i

        allocate(this%op_l,source=op_l)
        allocate(this%op_r,source=op_r)
        this%dim_mode=1
        do i=1,size(this%op_l)
            this%dim_mode(1)=this%dim_mode(1)*lat%get_order_dim(this%op_l(i))
        enddo
        do i=1,size(this%op_r)
            this%dim_mode(2)=this%dim_mode(2)*lat%get_order_dim(this%op_r(i))
        enddo
        this%set=.true.
    end subroutine

    subroutine init_otherH(this,H_in)
        class(t_H_base),intent(inout)    :: this
        class(t_H_base),intent(in)       :: h_in

        this%dimH=H_in%dimH
        this%op_l=H_in%op_l
        this%op_r=H_in%op_r
        this%dim_mode=H_in%dim_mode
        this%mult_M_single=H_in%mult_M_single
        this%desc=H_in%desc
        this%set=.true.
        this%dim_l_single=H_in%dim_l_single
        this%dim_r_single=H_in%dim_r_single
        if(allocated(H_in%mode_l)) Call H_in%mode_l%copy(this%mode_l)
        if(allocated(H_in%mode_r)) Call H_in%mode_r%copy(this%mode_r)
    end subroutine

    subroutine finish_setup_base(this)
        class(t_H_base),intent(inout)    :: this
        integer     ::  i

        if(.not.this%set) ERROR STOP "Hamiltoninan should be set when calling finish_setup"
        if(.not.allocated(this%mode_l)) ERROR STOP "mode_l of Hamiltoninan should be set when calling finish_setup"
        if(.not.allocated(this%mode_r)) ERROR STOP "mode_l of Hamiltoninan should be set when calling finish_setup"
        do i=1,number_different_order_parameters
            Call this%mode_l%get_mode_single_size(i,this%dim_l_single(i))
            Call this%mode_r%get_mode_single_size(i,this%dim_r_single(i))
        enddo
    end subroutine

    function is_set(this) result(l)
        class(t_H_base),intent(in)       ::  this
        logical                     ::  l

        l=this%set
    end function
    
    subroutine set_prepared(this,l)
        class(t_H_base),intent(inout)    ::  this
        logical,intent(in)          ::  l

        this%set=l
    end subroutine

    function same_space(this,comp) result(same)
        !checks if the Hamiltonian acts on the same space
        class(t_H_base),intent(in)    ::  this
        class(t_H_base),intent(in)    ::  comp
        logical                  ::  same

        same=.false.
        if(this%set.and.comp%set)then
            if(all(this%dimH==comp%dimH).and.all(this%dim_mode==comp%dim_mode))then
                if(size(this%op_l)==size(comp%op_l).and.size(this%op_r)==size(comp%op_r))then
                    if(all(this%op_l==comp%op_l).and.any(this%op_r==comp%op_r))then
                        if(this%mode_l%is_same(comp%mode_l))then
                            if(this%mode_r%is_same(comp%mode_r))then
                                same=.true.
                            endif
                        endif
                    endif
                endif
            endif
        endif
    end function



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!        ROUTINES THAT SHOULD BE IMPLEMENTED MORE EFFICIENTLY    !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_work_single(this,work,order)
    class(t_H_base),intent(inout)           :: this
    class(work_ham_single),intent(inout)    :: work 
    integer,intent(in)                      :: order

    ERROR STOP "IMPLEMENT"  !implement for local entries, then make deferred if I decide to keep that
end subroutine

subroutine set_work_mode(this,work)
    class(t_H_base),intent(inout)   :: this
    class(work_mode),intent(inout)  :: work 

    integer     :: sizes(2)

    sizes=0
    sizes(1)=this%dimH(1)+this%dimH(2)
    Call work%set(sizes)
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!            MPI ROUTINES           !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine send_base(this,ithread,tag,com)
    use mpi_basic                
    class(t_H_base),intent(in)      :: this
    integer,intent(in)              :: ithread
    integer,intent(in)              :: tag
    integer,intent(in)              :: com

#ifdef CPP_MPI
    integer     :: ierr
    integer     :: N(2),ident(2)

    if(.not.this%is_set()) ERROR STOP "CANNOT SEND HAMILTONIAN WHEN THE ORIGIN IS NOT SET"


    N=[size(this%op_l), size(this%op_r)]
    Call MPI_SEND(this%dimH,          size(this%dimH), MPI_INTEGER,   ithread, tag, com, ierr)
    Call MPI_SEND(N,                  2,               MPI_INTEGER,   ithread, tag, com, ierr)
    Call MPI_SEND(this%op_l,          N(1),            MPI_INTEGER,   ithread, tag, com, ierr)
    Call MPI_SEND(this%op_r,          N(2),            MPI_INTEGER,   ithread, tag, com, ierr)
    Call MPI_SEND(this%dim_mode,      2,               MPI_INTEGER,   ithread, tag, com, ierr)
    Call MPI_SEND(this%set,           1,               MPI_LOGICAL,   ithread, tag, com, ierr)
    Call MPI_SEND(this%desc,          len_desc,        MPI_CHARACTER, ithread, tag, com, ierr)
    Call MPI_SEND(this%mult_M_single, 1,               MPI_INTEGER,   ithread, tag, com, ierr)

    Call get_mode_ident(this%mode_l,ident(1))
    Call get_mode_ident(this%mode_r,ident(2))
    Call MPI_SEND(ident, 2, MPI_INTEGER, ithread, tag, com, ierr)

    Call this%mode_l%send(ithread,tag,com)
    Call this%mode_r%send(ithread,tag,com)

    Call MPI_SEND(this%dim_l_mode, number_different_order_parameters, MPI_INTEGER, ithread, tag, com, ierr)
    Call MPI_SEND(this%dim_r_mode, number_different_order_parameters, MPI_INTEGER, ithread, tag, com, ierr)
#else
    continue
#endif
end subroutine

subroutine recv_base(this,ithread,tag,com)
    use mpi_basic                
    class(t_H_base),intent(inout)   :: this
    integer,intent(in)              :: ithread
    integer,intent(in)              :: tag
    integer,intent(in)              :: com

#ifdef CPP_MPI
    integer     :: ierr
    integer     :: N(2),ident(2)
    integer     :: stat(MPI_STATUS_SIZE) 

    if(this%is_set()) ERROR STOP "CANNOT RECV HAMILTONIAN WHEN THE ORIGIN IS ALREADY SET"

    Call MPI_RECV(this%dimH,          size(this%dimH), MPI_INTEGER,   ithread, tag, com, stat, ierr)
    Call MPI_RECV(N,                  2,               MPI_INTEGER,   ithread, tag, com, stat, ierr)
    allocate(this%op_l(N(1)),this%op_r(N(2)))
    Call MPI_RECV(this%op_l,          N(1),            MPI_INTEGER,   ithread, tag, com, stat, ierr)
    Call MPI_RECV(this%op_r,          N(2),            MPI_INTEGER,   ithread, tag, com, stat, ierr)
    Call MPI_RECV(this%dim_mode,      2,               MPI_INTEGER,   ithread, tag, com, stat, ierr)
    Call MPI_RECV(this%set,           1,               MPI_LOGICAL,   ithread, tag, com, stat, ierr)
    Call MPI_RECV(this%desc,          len_desc,        MPI_CHARACTER, ithread, tag, com, stat, ierr)
    Call MPI_RECV(this%mult_M_single, 1,               MPI_INTEGER,   ithread, tag, com, stat, ierr)

    Call MPI_RECV(ident, 2, MPI_INTEGER, ithread, tag, com, stat, ierr)
    Call set_mode_ident(this%mode_l,ident(1))
    Call set_mode_ident(this%mode_r,ident(2))

    Call this%mode_l%recv(ithread,tag,com)
    Call this%mode_r%recv(ithread,tag,com)
    Call MPI_RECV(this%dim_l_single, number_different_order_parameters, MPI_INTEGER, ithread, tag, com, stat, ierr)
    Call MPI_RECV(this%dim_r_single, number_different_order_parameters, MPI_INTEGER, ithread, tag, com, stat, ierr)
    Call this%finish_setup()
#else
    continue
#endif
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!          HELPER ROUTINES          !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_mode_same(mode1,mode2)
    class(F_mode),intent(in),allocatable    ::  mode1,mode2
        if(allocated(mode1))then
            if(allocated(mode2))then
                if(.not.mode1%is_same(mode2)) ERROR STOP "Hamiltonian mode1 and mode2 differ"
            else
                ERROR STOP "mode1 allocated but mode2 not allocated, not sure if addition is in correct space"
            endif
        endif
end subroutine


subroutine get_sum(N1,N2,vec_in,vec_out)
    integer,intent(in)  ::  N1,N2
    real(8),intent(in)  :: vec_in(N1,N2)
    real(8),intent(out) :: vec_out(N2)
    
    vec_out=sum(vec_in,1)
end subroutine

subroutine mult_l_disc(this,i_m,lat,N,ind_out,vec,ind_sum,ind_Mult,mat_mult,vec_mult)
    !Calculates the entries of the left vector * matrix product for the indices ind_out of the result vector
    USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
    ! input
    class(t_H_base), intent(in)         :: this
    type(lattice), intent(in)           :: lat
    integer, intent(in)                 :: i_m          !index of the comp's right mode in the inner dim_mode
    integer, intent(in)                 :: N            !number of indices to calculated
    integer, intent(in)                 :: ind_out(N)   !indices to be calculated
    ! output
    real(8),intent(out)                 :: vec(N) ! dim_modes_inner(this%mode_r%order(comp))
    !temporary data
    integer,intent(inout)               :: ind_sum (N+1)                 !ind_mult index where the a vec-entry start and end
    integer,intent(inout)               :: ind_mult(N*this%col_max)     !indices of the left array which have non-vanishing contributions to get the ind entries of the vec/mat product
    real(8),intent(inout)               :: mat_mult(N*this%col_max)     !matrix entries corresponding to ind_mult
    real(8),intent(inout)               :: vec_mult(N*this%col_max)        !values of discontiguous left mode which has to be evaluated (indices of ind_mult)

    write(error_unit,'(//A)') "Trying to call mult_l_dist with Hamiltonian mode where this is not implemented."
    write(error_unit,'(A)')   "Please choose a Hamiltonian implementation including an implementation such as the eigen and mkl options."
    ERROR STOP
end subroutine 

subroutine mult_r_disc(this,i_m,lat,N,ind_out,vec,ind_sum,ind_Mult,mat_mult,vec_mult)
    !Calculates the entries of the matrix * right vector product for the indices ind_out of the result vector
    ! input
    class(t_H_base), intent(in)         :: this
    type(lattice), intent(in)           :: lat
    integer, intent(in)                 :: i_m          !index of the comp's right mode in the inner dim_mode
    integer, intent(in)                 :: N            !number of indices to calculated
    integer, intent(in)                 :: ind_out(N)   !indices to be calculated
    ! output
    real(8),intent(out)                 :: vec(N)
    !temporary data
    integer,intent(inout)               :: ind_sum (N+1)                !ind_mult index where the a vec-entry start and end
    integer,intent(inout)               :: ind_mult(N*this%row_max)     !indices of the right array which have non-vanishing contributions to get the ind entries of the mat/vecduct
    real(8),intent(inout)               :: mat_mult(N*this%row_max)     !matrix entries corresponding to ind_mult
    real(8),intent(inout)               :: vec_mult(N*this%row_max)     !values of discontiguous right mode which has to be evaluated (indices of ind_mult)

    write(error_unit,'(//A)') "Trying to call mult_l_dist with Hamiltonian mode where this is not implemented."
    write(error_unit,'(A)')   "Please choose a Hamiltonian implementation including an implementation such as the eigen and mkl options."
    ERROR STOP
end subroutine 
end module
