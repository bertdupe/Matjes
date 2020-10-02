module m_H_type
!module containing the basic polymorphic Hamiltonian class t_H
use m_derived_types, only : lattice
use m_derived_types, only: operator_real_order_N
implicit none

type,abstract :: t_H
    integer,allocatable     :: op_l(:),op_r(:)
    logical,private         :: set=.false.
contains
    procedure  :: is_set
    procedure  :: set_prepared
    procedure(int_eval_single),deferred :: eval_single
    procedure(int_eval_all),deferred    :: eval_all
    procedure(int_set_H),deferred       :: set_H
    procedure(int_set_H_1),deferred     :: set_H_1
    procedure(int_add_H),deferred       :: add_H
    procedure(int_destroy),deferred     :: destroy
    procedure(int_copy),deferred        :: copy
    procedure(int_destroy),deferred     :: optimize
    procedure(int_mult),deferred        :: mult_r,mult_l


    procedure,NON_OVERRIDABLE           :: destroy_base
    procedure,NON_OVERRIDABLE           :: add_base
    procedure,NON_OVERRIDABLE           :: copy_base
end type
private
public t_H,energy_all

interface
    subroutine int_mult(this,lat,vec)
        import t_H,lattice
    	class(t_H),intent(in)     :: this
    	type(lattice),intent(in)  :: lat
        real(8),intent(inout)     :: vec(:)
    end subroutine

    subroutine int_destroy(this)
        import t_H
    	class(t_H),intent(inout)  :: this
    end subroutine

    subroutine int_set_H(this,energy_in,lat)
        import t_H,lattice, operator_real_order_N
    	!need to get rid of dim_mode input
    	class(t_H),intent(inout)  :: this
    	type(operator_real_order_N)     :: energy_in
    	type(lattice),intent(in)    	:: lat
    end subroutine

    subroutine int_eval_all(this,E, lat)
        import t_H,lattice
        class(t_H),intent(in)    ::  this
        type(lattice),intent(in)    ::  lat
        real(8),intent(out)			::	E
    end subroutine

    subroutine int_eval_single(this,E,i_m, lat)
        import t_H,lattice
        class(t_H),intent(in)    ::  this
        type(lattice),intent(in)    ::  lat
    	integer,intent(in)			::	i_m
        real(8),intent(out)			::	E
    end subroutine
    subroutine int_set_H_1(this,line,Hval,Hval_ind,order,lat)
        import t_h,lattice
        class(t_H),intent(inout)    :: this
    	type(lattice),intent(in)    :: lat
        integer,intent(in)          :: order(2)
        real(8),intent(in)          :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
        integer,intent(in)          :: Hval_ind(:,:)
        integer,intent(in)          :: line(:,:)
    end subroutine

    subroutine int_add_H(this,H_add)
        import t_h
        class(t_H),intent(inout)    :: this
        class(t_H),intent(in)       :: H_add
    end subroutine

    subroutine int_copy(this,Hout)
        import t_h
        class(t_H),intent(in)         :: this
        class(t_H),intent(inout)      :: Hout
    end subroutine

end interface

contains


    subroutine copy_base(this,Hout)
        class(t_H),intent(in)         :: this
        class(t_H),intent(inout)      :: Hout

        if(this%is_set())then
            if(.not.allocated(this%op_l).or..not.allocated(this%op_r))&
                & STOP "Cannot copy H since base op_l/op_r are not allocated"
            allocate(Hout%op_l,source=this%op_l)
            allocate(Hout%op_r,source=this%op_r)
            Call Hout%set_prepared(.true.)
        else
            STOP "cannot copy H since source is not set"
        endif
    end subroutine

    subroutine add_base(this,add)
        class(t_H),intent(inout)    :: this
        class(t_H),intent(in)       :: add

        if(add%is_set())then
            if(.not.allocated(add%op_l).or..not.allocated(add%op_r)) &
                & STOP "cannot add H and its op_l/op_r-arguments are not allocated"
            allocate(this%op_l,source=add%op_l)
            allocate(this%op_r,source=add%op_r)
            Call this%set_prepared(.true.)
        else
            STOP "cannot add H since argument is not set"
        endif
    end subroutine

    subroutine destroy_base(this)
        class(t_H),intent(inout)    :: this

        if(allocated(this%op_l)) deallocate(this%op_l)
        if(allocated(this%op_r)) deallocate(this%op_r)
        Call this%set_prepared(.false.)
    end subroutine

    subroutine energy_all(ham,lat,E)
        class(t_H),intent(in)       ::  ham(:)
        class(lattice),intent(in)   ::  lat
        real(8),intent(out)         ::  E

        real(8)     ::  tmp_E(size(ham))
        integer     ::  i

        E=0.0d0
        do i=1,size(ham)
          Call ham(i)%eval_all(tmp_E(i),lat)
        enddo
        E=sum(tmp_E)
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


end module
