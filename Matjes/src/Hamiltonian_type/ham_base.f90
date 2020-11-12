module m_H_type
!module containing the basic polymorphic Hamiltonian class t_H
use m_derived_types, only : lattice
use m_derived_types, only: operator_real_order_N
implicit none

type,abstract :: t_H
    integer                 :: dimH(2)=0 !dimension of Hamiltonian
    integer,allocatable     :: op_l(:),op_r(:)  !operator indices which have to be combined to get the space of the left/right side of the Hamiltonian
    integer                 :: dim_mode(2)=0     !size of left/right mode after multiplying out order parameters (size per lattice point)
    logical,private         :: set=.false. !has this object been set?
contains

    !Hamiltonian initialization routines
    procedure(int_init_H),deferred           :: init !based on old operator_real_order_N input
    procedure(int_init_H_1),deferred         :: init_1 !based on Hamiltonian in coo format input (arrays), only rank_2 
    procedure(int_init_H_mult_2),deferred    :: init_mult_2 !for rank >2, but only at 2 sites at once

    procedure(int_destroy),deferred         :: optimize  !calls possible optimization routines rearanging internal Hamiltonian layout
    procedure(int_mult),deferred            :: mult_r,mult_l !multipy out with left/right side

    !evaluation routines
        !actually these could be overridable, but so far used to check no old implementation still has these wrongly
    procedure ,NON_OVERRIDABLE              :: eval_all !evaluates energy of full lattice
    procedure ,NON_OVERRIDABLE              :: mult_r_red,mult_l_red !multiplied out left/right and reduces result to only one order-parameter
    procedure ,NON_OVERRIDABLE              :: mult_r_red_single,mult_l_red_single

        !should be implemented more efficiently, where possible
    procedure                               :: mult_r_single,mult_l_single !multipy out with left/right side

    !really non_overridable
    procedure,NON_OVERRIDABLE               :: destroy
    procedure,NON_OVERRIDABLE               :: copy
    procedure,NON_OVERRIDABLE               :: add
    procedure,NON_OVERRIDABLE               :: init_base

    procedure(int_add_H),deferred           :: add_child
    procedure(int_destroy),deferred         :: destroy_child
    procedure(int_copy),deferred            :: copy_child

    !misc.
    procedure,NON_OVERRIDABLE  :: is_set
    procedure,NON_OVERRIDABLE  :: set_prepared

    !TODO:
    procedure(int_eval_single),deferred     :: eval_single !needs some work

    !finalize routine? (might be risky with Hamiltonian references that have been passed around)
end type
private
public t_H

interface
    subroutine int_mult(this,lat,res)
        import t_H,lattice
        class(t_H),intent(in)     :: this
        type(lattice),intent(in)  :: lat
        real(8),intent(inout)     :: res(:)
    end subroutine


    subroutine int_mult_single(this,i_site,lat,res)
        import t_H,lattice
        class(t_H),intent(in)     :: this
        integer,intent(in)        :: i_site
        type(lattice),intent(in)  :: lat
        real(8),intent(inout)     :: res(:)
    end subroutine

    subroutine int_mult_red(this,lat,res,op_keep)
        import t_H,lattice
        class(t_H),intent(in)     :: this
        type(lattice),intent(in)  :: lat
        real(8),intent(inout)     :: res(:)
        integer,intent(in)        :: op_keep
    end subroutine

    subroutine int_destroy(this)
        import t_H
        class(t_H),intent(inout)  :: this
    end subroutine

    subroutine int_init_H(this,energy_in,lat)
        import t_H,lattice, operator_real_order_N
        !need to get rid of dim_mode input
        class(t_H),intent(inout)  :: this
        type(operator_real_order_N)     :: energy_in
        type(lattice),intent(in)        :: lat
    end subroutine

    subroutine int_eval_all(this,E, lat)
        import t_H,lattice
        class(t_H),intent(in)    ::  this
        type(lattice),intent(in)    ::  lat
        real(8),intent(out)         ::  E
    end subroutine

    subroutine int_eval_single(this,E,i_m, lat)
        import t_H,lattice
        class(t_H),intent(in)    ::  this
        type(lattice),intent(in)    ::  lat
        integer,intent(in)          ::  i_m
        real(8),intent(out)         ::  E
    end subroutine
    subroutine int_init_H_1(this,line,Hval,Hval_ind,order,lat)
        import t_h,lattice
        class(t_H),intent(inout)    :: this
        type(lattice),intent(in)    :: lat
        integer,intent(in)          :: order(2)
        real(8),intent(in)          :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
        integer,intent(in)          :: Hval_ind(:,:)
        integer,intent(in)          :: line(:,:)
    end subroutine

    subroutine int_init_H_mult_2(this,connect,Hval,Hval_ind,op_l,op_r,lat)
        !Constructs a Hamiltonian that depends on more than 2 order parameters but only at 2 sites (i.e. some terms are onsite)
        !(example: ME-coupling M_i*E_i*M_j
        import t_H,lattice
        class(t_H),intent(inout)    :: this
        type(lattice),intent(in)        :: lat
        !input Hamiltonian
        real(8),intent(in)              :: Hval(:)  !values of local Hamiltonian for each line
        integer,intent(in)              :: Hval_ind(:,:)  !indices in order-parameter space for Hval
        integer,intent(in)              :: op_l(:),op_r(:) !which order parameters are used at left/right side of local Hamiltonian-matrix
        integer,intent(in)              :: connect(:,:) !lattice sites to be connected (2,Nconnections)
    end subroutine


    subroutine int_add_H(this,H_in)
        import t_h
        class(t_H),intent(inout)    :: this
        class(t_H),intent(in)       :: H_in
    end subroutine

    subroutine int_copy(this,Hout)
        import t_h
        class(t_H),intent(in)         :: this
        class(t_H),intent(inout)      :: Hout
    end subroutine

end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!                    MOST USEFULL ROUTINES                        !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine eval_all(this,E,lat)
        USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_DOUBLE
        !evaluates the energie summed for all sites
        class(t_H),intent(in)       :: this
        type(lattice), intent(in)   :: lat
        real(8), intent(out)        :: E
        ! internal
        real(C_DOUBLE)              :: tmp(this%dimH(1))
        real(8),pointer             :: modes_l(:)
        real(8),allocatable,target  :: vec_l(:)
#ifdef CPP_BLAS
        real(8),external            :: ddot !blas routine
#endif
    
        Call lat%point_order(this%op_l,this%dimH(1),modes_l,vec_l)
    
        Call this%mult_r(lat,tmp)
#ifdef CPP_BLAS
        E=ddot(this%dimH(1),modes_l,1,tmp,1)
#else
        E=dot_product(modes_l,tmp)
#endif
    
        if(allocated(vec_l)) deallocate(vec_l) 
    end subroutine 



    subroutine mult_r_red(this,lat,res,op_keep)
        !multiply out right side, multiply with left order parameter, reduce to only keep operator corresponding to op_keep
        !this is mainly necessary to calculate the effective magnetic field (corresponding to derivative with respect to one order parameter)
        use m_derived_types, only: lattice
        class(t_H),intent(in)           :: this
        type(lattice), intent(in)       :: lat
        real(8), intent(inout)          :: res(:)   !result matrix-vector product
        integer,intent(in)              :: op_keep
        ! internal
        real(8),allocatable             :: tmp(:)   !multipied, but not reduced
        real(8),pointer                 :: modes(:)
        real(8),allocatable,target      :: vec(:)
    
        allocate(tmp(this%dimH(1)))
        Call this%mult_r(lat,tmp)
        allocate(vec(this%dimH(1)),source=0.0d0)
        Call lat%set_order_comb_exc(this%op_l,vec,this%op_l==op_keep)
        tmp=tmp*vec
        Call lat%reduce(tmp,this%op_l,op_keep,res)
        deallocate(tmp,vec)
    end subroutine 
    
    subroutine mult_l_red(this,lat,res,op_keep)
        !multiply out left side, multiply with right order parameter, reduce to only keep operator corresponding to op_keep
        !this is mainly necessary to calculate the effective magnetic field (corresponding to derivative with respect to one order parameter)
        use m_derived_types, only: lattice
        class(t_H),intent(in)           :: this
        type(lattice), intent(in)       :: lat
        real(8), intent(inout)          :: res(:)   !result matrix-vector product
        integer,intent(in)              :: op_keep
        ! internal
        real(8),allocatable             :: tmp(:)   !multipied, but not reduced
        real(8),allocatable,target      :: vec(:)
    
        allocate(tmp(this%dimH(2)))
        Call this%mult_l(lat,tmp)
        allocate(vec(this%dimH(2)),source=0.0d0)
        Call lat%set_order_comb_exc(this%op_r,vec,this%op_r==op_keep)
        tmp=tmp*vec
        Call lat%reduce(tmp,this%op_r,op_keep,res)
        deallocate(tmp,vec)
    end subroutine 

    subroutine mult_r_red_single(this,i_site,lat,res,op_keep)
        !multiply out right side, multiply with left order parameter, reduce to only keep operator corresponding to op_keep
        !this is mainly necessary to calculate the effective magnetic field (corresponding to derivative with respect to one order parameter)
        use m_derived_types, only: lattice
        class(t_H),intent(in)           :: this
        type(lattice), intent(in)       :: lat
        integer,intent(in)              :: i_site
        real(8), intent(inout)          :: res(:)   !result matrix-vector product
        integer,intent(in)              :: op_keep
        ! internal
        real(8),allocatable             :: tmp(:)   !multipied, but not reduced
        real(8),pointer                 :: modes(:)
        real(8),allocatable,target      :: vec(:)
    
        allocate(tmp(this%dim_mode(1)))
        Call this%mult_r_single(i_site,lat,tmp)
        allocate(vec(this%dim_mode(1)),source=0.0d0)
        Call lat%set_order_comb_exc_single(i_site,this%op_l,vec,this%op_l==op_keep)
        tmp=tmp*vec
        Call lat%reduce_single(i_site,tmp,this%op_l,op_keep,res)
        deallocate(tmp,vec)
    end subroutine 
    
    subroutine mult_l_red_single(this,i_site,lat,res,op_keep)
        !multiply out left side, multiply with right order parameter, reduce to only keep operator corresponding to op_keep
        !this is mainly necessary to calculate the effective magnetic field (corresponding to derivative with respect to one order parameter)
        use m_derived_types, only: lattice
        class(t_H),intent(in)           :: this
        integer,intent(in)              :: i_site
        type(lattice), intent(in)       :: lat
        real(8), intent(inout)          :: res(:)   !result matrix-vector product
        integer,intent(in)              :: op_keep
        ! internal
        real(8),allocatable             :: tmp(:)   !multipied, but not reduced
        real(8),allocatable,target      :: vec(:)
    
        allocate(tmp(this%dim_mode(2)))
        Call this%mult_l_single(i_site,lat,tmp)
        allocate(vec(this%dim_mode(2)),source=0.0d0)
        Call lat%set_order_comb_exc_single(i_site,this%op_r,vec,this%op_r==op_keep)
        tmp=tmp*vec
        Call lat%reduce_single(i_site,tmp,this%op_r,op_keep,res)
        deallocate(tmp,vec)
    end subroutine 



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!                    BASE ROUTINES                                !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine copy(this,Hout)
        !copy this to Hout
        class(t_H),intent(in)         :: this
        class(t_H),intent(inout)      :: Hout

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
        else
            STOP "cannot copy H since source is not set"
        endif
    end subroutine

    subroutine add(this,H_in)
        !add H_in to this Hamiltonian, or set this to H_in is this is not set
        class(t_H),intent(inout)    :: this
        class(t_H),intent(in)       :: H_in

        if(this%is_set())then
            if(.not.allocated(H_in%op_l).or..not.allocated(H_in%op_r)) &
                & STOP "cannot add H and its op_l/op_r-arguments are not allocated"
            if(.not.all(this%op_l==H_in%op_l)) STOP "CANNOT ADD hamiltonians with different op_l"
            if(.not.all(this%op_r==H_in%op_r)) STOP "CANNOT ADD hamiltonians with different op_r"
            if(any(this%dimH/=H_in%dimH)) STOP "Trying to add Hamiltonians with different Hamiltonian dimensions"
            Call this%add_child(H_in)
        else
            Call H_in%copy(this)
        endif
    end subroutine

    subroutine destroy(this)
        !Destroys Hamiltonian
        class(t_H),intent(inout)    :: this

        Call this%destroy_child()
        this%dimH=0
        if(allocated(this%op_l)) deallocate(this%op_l)
        if(allocated(this%op_r)) deallocate(this%op_r)
        this%dim_mode=0
        Call this%set_prepared(.false.)
    end subroutine

    subroutine init_base(this,lat,op_l,op_r)
        class(t_H),intent(inout)    :: this
        class(lattice),intent(in)   :: lat
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

    function is_set(this) result(l)
        class(t_H),intent(in)       ::  this
        logical                     ::  l

        l=this%set
    end function
    
    subroutine set_prepared(this,l)
        class(t_H),intent(inout)    ::  this
        logical,intent(in)          ::  l

        this%set=l
    end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!        ROUTINES THAT SHOULD BE IMPLEMENTED MORE EFFICIENTLY    !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mult_l_single(this,i_site,lat,res)
    use m_derived_types, only: lattice
    class(t_H),intent(in)        :: this
    integer,intent(in)           :: i_site
    type(lattice),intent(in)     :: lat
    real(8),intent(inout)        :: res(:)

    real(8)                      :: tmp(this%dimH(2))

    Call this%mult_l(lat,tmp)
    res=tmp((i_site-1)*this%dim_mode(1)+1:i_site*this%dim_mode(1))

end subroutine

subroutine mult_r_single(this,i_site,lat,res)
    use m_derived_types, only: lattice
    class(t_H),intent(in)        :: this
    integer,intent(in)           :: i_site
    type(lattice),intent(in)     :: lat
    real(8),intent(inout)        :: res(:)

    real(8)                      :: tmp(this%dimH(1))

    Call this%mult_r(lat,tmp)
    res=tmp((i_site-1)*this%dim_mode(2)+1:i_site*this%dim_mode(2))
end subroutine


end module
