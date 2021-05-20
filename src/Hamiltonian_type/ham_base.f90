module m_H_type
!module containing the basic polymorphic Hamiltonian class t_H_base
use m_derived_types, only : lattice,number_different_order_parameters
use m_mode_public, only: F_mode, get_mode_ident, set_mode_ident
use m_work_ham_single, only:  work_ham_single
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
contains

    !Hamiltonian initialization routines
    procedure(int_init_H_connect),deferred          :: init_connect         !based on Hamiltonian in coo format input (arrays), only rank_2, connection input
    procedure(int_init_H_mult_connect_2),deferred   :: init_mult_connect_2  !based on Hamiltonian in coo format input (arrays), for rank >2, but only at 2 sites at once, connection input
    procedure(int_init_H_coo),deferred              :: init_coo             !based on Hamiltonian in coo format input (arrays), the coo data-arrays are directly moved into the type
    procedure(int_destroy),deferred                 :: optimize             !calls possible optimization routines rearanging internal Hamiltonian layout

    !routines acting on the entire Hamiltonian
    procedure ,NON_OVERRIDABLE              :: eval_all                 !evaluates energy of full lattice
    procedure(int_mult),deferred            :: mult_r,mult_l            !multipy out with left/right side
    procedure ,NON_OVERRIDABLE              :: mult_r_red,mult_l_red    !multiplied out left/right and reduces result to only one order-parameter
    procedure ,NON_OVERRIDABLE              :: energy_dist              !get energy distribution per site

    !routines acting on parts of the Hamiltonian
    procedure                               :: mult_l_ind,      mult_r_ind          !
    procedure                               :: mult_l_disc_disc,mult_r_disc_disc
    procedure                               :: get_ind_mult_r, get_ind_mult_l       !get the indices that appear mulitplying a sparse vectors with the matrix
!    !probably never used anywhere,implemented for eigen, keep so far as it might make sense at some point
!    procedure(int_mult_cont),deferred       :: mult_l_cont,mult_r_cont  
!    procedure(int_mult_disc),deferred       :: mult_l_disc,mult_r_disc

    !routines acting on parts of the Hamiltonian defined by a single site of an order-parameter
    procedure(int_eval_single),deferred     :: eval_single                              !evaluates the energy of a single site of a given order-parameter
    procedure                               :: eval_single_work                         !evaluates the energy caused by a single mode (using a work array)
    procedure                               :: set_work_size_single                     !sets the necessary sizes for the work arrays

    procedure                               :: mult_l_single,   mult_r_single           !multipy out with left/right side
    procedure ,NON_OVERRIDABLE              :: mult_r_red_single,mult_l_red_single      !multipy out with left/right side excluding single order-parameter

    !Utility functions
    procedure,NON_OVERRIDABLE               :: destroy
    procedure,NON_OVERRIDABLE               :: copy
    procedure,NON_OVERRIDABLE               :: add
    procedure,NON_OVERRIDABLE               :: init_base
    procedure,NON_OVERRIDABLE               :: init_otherH
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
    subroutine int_mult(this,lat,res)
        import t_H_base,lattice
        class(t_H_base),intent(in)     :: this
        type(lattice),intent(in)  :: lat
        real(8),intent(inout)     :: res(:)
    end subroutine

    subroutine int_mult_single(this,i_site,lat,res)
        import t_H_base,lattice
        class(t_H_base),intent(in)     :: this
        integer,intent(in)        :: i_site
        type(lattice),intent(in)  :: lat
        real(8),intent(inout)     :: res(:)
    end subroutine

    subroutine int_mult_red(this,lat,res,op_keep)
        import t_H_base,lattice
        class(t_H_base),intent(in)     :: this
        type(lattice),intent(in)  :: lat
        real(8),intent(inout)     :: res(:)
        integer,intent(in)        :: op_keep
    end subroutine

    subroutine int_mult_cont(this,bnd,vec,res)
        !multiply to the either side with a continuous section of the applied vector
        import t_H_base
        class(t_H_base),intent(in)       :: this
        integer,intent(in)          :: bnd(2)
        real(8),intent(in)          :: vec(bnd(2)-bnd(1)+1)
        real(8),intent(inout)       :: res(:)   !result matrix-vector product
    end subroutine 
    
    subroutine int_mult_disc(this,N,ind,vec,res)
        !multiply to the either side with a discontinuous section of the applied vector
        import t_H_base
        class(t_H_base),intent(in)       :: this
        integer,intent(in)          :: N
        integer,intent(in)          :: ind(N)
        real(8),intent(in)          :: vec(N)
        real(8),intent(inout)       :: res(:)   !result matrix-vector product
    end subroutine 

    subroutine int_mult_ind(this,vec,N,ind_out,vec_out)
        import t_H_base
        class(t_H_base),intent(in)       :: this
        real(8),intent(in)          :: vec(:)
        integer,intent(in)          :: N
        integer,intent(in)          :: ind_out(N)
        real(8),intent(out)         :: vec_out(N)
    end subroutine 


    subroutine int_destroy(this)
        import t_H_base
        class(t_H_base),intent(inout)  :: this
    end subroutine

    subroutine int_eval_all(this,E, lat)
        import t_H_base,lattice
        class(t_H_base),intent(in)    ::  this
        type(lattice),intent(in)    ::  lat
        real(8),intent(out)         ::  E
    end subroutine

    subroutine int_eval_single(this,E,i_m,order,lat)
        import t_H_base,lattice,number_different_order_parameters
        class(t_H_base),intent(in)    ::  this
        type(lattice),intent(in)    ::  lat
        integer,intent(in)          ::  i_m     !site index in (1:Ncell*dim_mode_inner(order)/dim_mode(order))-basis
        integer,intent(in)          :: order    
!        integer,intent(in)          ::  dim_bnd(2,number_different_order_parameters)    !starting/final index in respective dim_mode of the order parameter (so that energy of single magnetic atom can be be calculated
        real(8),intent(out)         ::  E       !energy caused by considered site
    end subroutine

    subroutine int_init_H_mult_connect_2(this,connect,Hval,Hval_ind,op_l,op_r,lat,mult_M_single,dim_mode_in)
        import t_H_base,lattice
        class(t_H_base),intent(inout)    :: this
        type(lattice),intent(in)    :: lat
        character(*),intent(in)     :: op_l
        character(*),intent(in)     :: op_r
        real(8),intent(in)          :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
        integer,intent(in)          :: Hval_ind(:,:)
        integer,intent(in)          :: connect(:,:)
        integer,intent(in)          :: mult_M_single
        integer,intent(in),optional :: dim_mode_in(2)   !optional way of putting in dim_mode directly (mainly for custom(not fully unfolded)rankN tensors)
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
        class(t_H_base),intent(inout)    :: this
        type(lattice),intent(in)    :: lat
        character(2),intent(in)     :: order
        real(8),intent(in)          :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
        integer,intent(in)          :: Hval_ind(:,:)
        integer,intent(in)          :: connect(:,:)
        integer,intent(in)          :: mult_M_single
    end subroutine

    subroutine int_add_H(this,H_in)
        import t_H_base
        class(t_H_base),intent(inout)    :: this
        class(t_H_base),intent(in)       :: H_in
    end subroutine

    subroutine int_copy(this,Hout)
        import t_H_base
        class(t_H_base),intent(in)         :: this
        class(t_H_base),intent(inout)      :: Hout
    end subroutine

    subroutine int_mpicom(this,comm)
        use mpi_basic                
        import t_H_base
        class(t_H_base),intent(inout)   ::  this
        type(mpi_type),intent(in)       ::  comm
    end subroutine

    subroutine int_send(this,ithread,tag,com)
        import t_H_base
        class(t_H_base),intent(in)  :: this
        integer,intent(in)          :: ithread
        integer,intent(in)          :: tag
        integer,intent(in)          :: com
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

    subroutine eval_all(this,E,lat)
        USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_DOUBLE
        use, intrinsic :: iso_fortran_env, only : error_unit
        !evaluates the energie summed for all sites
        class(t_H_base),intent(in)       :: this
        type(lattice), intent(in)   :: lat
        real(8), intent(out)        :: E
        ! internal
        real(C_DOUBLE)              :: tmp(this%dimH(1))
        real(8),pointer             :: modes_l(:)
        real(8),allocatable,target  :: vec_l(:)
    
        if(.not.allocated(this%mode_l).or..not.allocated(this%mode_r))then
            write(error_unit,'(2/2A)') "Failed to evaluate Hamiltonian: ", this%desc
            ERROR STOP "IMPLEMENT mode_l/mode_r for all Hamiltonians"
        endif

        Call this%mode_l%get_mode(lat,modes_l,vec_l)
    
        Call this%mult_r(lat,tmp)
        E=dot_product(modes_l,tmp)

        if(allocated(vec_l)) deallocate(vec_l) 
    end subroutine 

    subroutine energy_dist(this,lat,order,E)
        use m_type_lattice, only: dim_modes_inner
        !get the energy values per site for a given order
        !might be wrong for non-symmetrical Hamiltonians
        class(t_H_base),intent(in)       :: this
        type(lattice), intent(in)   :: lat
        integer,intent(in)          :: order
        real(8),intent(out)         :: E(lat%Ncell*lat%site_per_cell(order))
        ! internal
        real(8),pointer             :: modes(:)
        real(8),allocatable,target  :: vec(:)
        real(8),allocatable         :: tmp(:)
        real(8),allocatable         :: red(:)

        if(any(this%op_l==order))then
            allocate(tmp(this%dimH(1)))
            Call this%mult_r(lat,tmp)
            Call this%mode_l%get_mode(lat,modes,vec)
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
            Call this%mult_l(lat,tmp)
            Call this%mode_r%get_mode(lat,modes,vec)
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
        if(allocated(vec)) deallocate(vec)
    end subroutine

    subroutine mult_r_red(this,lat,res,op_keep)
        !multiply out right side, multiply with left order parameter, reduce to only keep operator corresponding to op_keep
        !this is mainly necessary to calculate the effective magnetic field (corresponding to derivative with respect to one order parameter)
        class(t_H_base),intent(in)           :: this
        type(lattice), intent(in)       :: lat
        real(8), intent(inout)          :: res(:)   !result matrix-vector product
        integer,intent(in)              :: op_keep
        ! internal
        real(8)                         :: tmp(this%dimH(1))   !multipied, but not reduced
    
        Call this%mult_r(lat,tmp)
        Call this%mode_l%reduce_other_exc(lat,op_keep,tmp,res)
    end subroutine 
    
    subroutine mult_l_red(this,lat,res,op_keep)
        !multiply out left side, multiply with right order parameter, reduce to only keep operator corresponding to op_keep
        !this is mainly necessary to calculate the effective magnetic field (corresponding to derivative with respect to one order parameter)
        class(t_H_base),intent(in)           :: this
        type(lattice), intent(in)       :: lat
        real(8), intent(inout)          :: res(:)   !result matrix-vector product
        integer,intent(in)              :: op_keep
        ! internal
        real(8)                         :: tmp(this%dimH(2))   !multipied, but not reduced
    
        Call this%mult_l(lat,tmp)
        Call this%mode_r%reduce_other_exc(lat,op_keep,tmp,res)
    end subroutine 

    subroutine mult_l_ind(this,lat,N,ind_out,vec_out)
        class(t_H_base),intent(in)       :: this
        type(lattice),intent(in)    :: lat
        integer,intent(in)          :: N
        integer,intent(in)          :: ind_out(N)
        real(8),intent(out)         :: vec_out(N)

        real(8)                     :: res(this%dimH(2))

        Call this%mult_l(lat,res)
        vec_out=res(ind_out)
    end subroutine

    subroutine mult_r_ind(this,lat,N,ind_out,vec_out)
        class(t_H_base),intent(in)  :: this
        type(lattice),intent(in)    :: lat
        integer,intent(in)          :: N
        integer,intent(in)          :: ind_out(N)
        real(8),intent(out)         :: vec_out(N)

        real(8)                     :: res(this%dimH(1))

        Call this%mult_r(lat,res)
        vec_out=res(ind_out)
    end subroutine

    subroutine mult_r_red_single(this,i_site,lat,res,op_keep)
        !multiply out right side, multiply with left order parameter, reduce to only keep operator corresponding to op_keep
        !this is mainly necessary to calculate the effective magnetic field (corresponding to derivative with respect to one order parameter)
        class(t_H_base),intent(in)           :: this
        type(lattice), intent(in)       :: lat
        integer,intent(in)              :: i_site
        real(8), intent(inout)          :: res(:)   !result matrix-vector product
        integer,intent(in)              :: op_keep
        ! internal
        real(8),allocatable             :: tmp(:)   !multipied, but not reduced
        real(8),pointer                 :: modes(:)
        real(8),allocatable,target      :: vec(:)

        ERROR STOP "NOT UPDATED TO NEW MODE TYPE"
    
        allocate(tmp(this%dim_mode(1)))
        Call this%mult_r_single(i_site,lat,tmp)
        allocate(vec(this%dim_mode(1)),source=0.0d0)
        Call lat%set_order_comb_exc_single(i_site,this%op_l,vec,this%op_l==op_keep)
        tmp=tmp*vec
!        Call lat%reduce_single(i_site,tmp,this%op_l,op_keep,res)   !op_keep changed to ind_keep
        deallocate(tmp,vec)
    end subroutine 
    
    subroutine mult_l_red_single(this,i_site,lat,res,op_keep)
        !multiply out left side, multiply with right order parameter, reduce to only keep operator corresponding to op_keep
        !this is mainly necessary to calculate the effective magnetic field (corresponding to derivative with respect to one order parameter)
        class(t_H_base),intent(in)           :: this
        integer,intent(in)              :: i_site
        type(lattice), intent(in)       :: lat
        real(8), intent(inout)          :: res(:)   !result matrix-vector product
        integer,intent(in)              :: op_keep
        ! internal
        real(8),allocatable             :: tmp(:)   !multipied, but not reduced
        real(8),allocatable,target      :: vec(:)
    
        ERROR STOP "NOT UPDATED TO NEW MODE TYPE"

        allocate(tmp(this%dim_mode(2)))
        Call this%mult_l_single(i_site,lat,tmp)
        allocate(vec(this%dim_mode(2)),source=0.0d0)
        Call lat%set_order_comb_exc_single(i_site,this%op_r,vec,this%op_r==op_keep)
        tmp=tmp*vec
!        Call lat%reduce_single(i_site,tmp,this%op_r,op_keep,res)   !op_keep changed to ind_keep
        deallocate(tmp,vec)
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
        if(allocated(H_in%mode_l)) Call H_in%mode_l%copy(this%mode_l)
        if(allocated(H_in%mode_r)) Call H_in%mode_r%copy(this%mode_r)
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
subroutine mult_l_single(this,i_site,lat,res)
    class(t_H_base),intent(in)        :: this
    integer,intent(in)           :: i_site
    type(lattice),intent(in)     :: lat
    real(8),intent(inout)        :: res(:)

    real(8)                      :: tmp(this%dimH(2))

    Call this%mult_l(lat,tmp)
    res=tmp((i_site-1)*this%dim_mode(1)+1:i_site*this%dim_mode(1))
end subroutine

subroutine mult_r_single(this,i_site,lat,res)
    class(t_H_base),intent(in)   :: this
    integer,intent(in)           :: i_site
    type(lattice),intent(in)     :: lat
    real(8),intent(inout)        :: res(:)

    real(8)                      :: tmp(this%dimH(1))

    Call this%mult_r(lat,tmp)
    res=tmp((i_site-1)*this%dim_mode(2)+1:i_site*this%dim_mode(2))
end subroutine

subroutine get_ind_mult_r(this,ind_in,N_out,ind_out)
    class(t_H_base),intent(in)      :: this
    integer,intent(in)              :: ind_in(:)
    integer,intent(out)             :: N_out
    integer,intent(inout)           :: ind_out(:)

    ERROR STOP "IMPLEMENT LOCALLY"
end subroutine

subroutine get_ind_mult_l(this,ind_in,N_out,ind_out)
    class(t_H_base),intent(in)      :: this
    integer,intent(in)              :: ind_in(:)
    integer,intent(out)             :: N_out
    integer,intent(inout)           :: ind_out(:)

    ERROR STOP "IMPLEMENT LOCALLY"
end subroutine

subroutine mult_r_disc_disc(this,ind_r,vec_r,ind_l,vec_out)
    class(t_H_base),intent(in)      :: this
    integer,intent(in)              :: ind_r(:)
    real(8),intent(in)              :: vec_r(:)
    integer,intent(in)              :: ind_l(:)
    real(8),intent(inout)           :: vec_out(:)

    ERROR STOP "IMPLEMENT LOCALLY"
end subroutine

subroutine mult_l_disc_disc(this,ind_l,vec_l,ind_r,vec_out)
    class(t_H_base),intent(in)           :: this
    integer,intent(in)              :: ind_l(:)
    real(8),intent(in)              :: vec_l(:)
    integer,intent(in)              :: ind_r(:)
    real(8),intent(inout)           :: vec_out(:)

    ERROR STOP "IMPLEMENT LOCALLY"
end subroutine

subroutine eval_single_work(this,E,i_m,order,lat,work)
    use m_derived_types, only: lattice, dim_modes_inner
    USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
    ! input
    class(t_H_base), intent(in)         :: this
    type(lattice), intent(in)           :: lat
    integer, intent(in)                 :: i_m
    integer, intent(in)                 :: order
    ! output
    real(8), intent(out)                :: E
    !temporary data
    type(work_ham_single),intent(inout) ::  work

    ERROR STOP "IMPLEMENT"  !implement for local entries, then make deferred if I decide to keep that
end subroutine

subroutine set_work_size_single(this,work,order)
    class(t_H_base),intent(inout)           :: this
    class(work_ham_single),intent(inout)    :: work 
    integer,intent(in)                      :: order

    ERROR STOP "IMPLEMENT"  !implement for local entries, then make deferred if I decide to keep that
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
end module
