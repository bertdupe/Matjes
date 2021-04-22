module m_H_coo_based
use m_H_deriv, only: t_H
use m_H_type, only: t_H_base
use m_derived_types, only: lattice
use m_H_coo, only: t_H_coo

public

type,abstract,extends(t_H) :: t_H_coo_based
    contains
    procedure(int_set_from_Hcoo),deferred   :: set_from_Hcoo    !routine to set the local parameters of the Hamiltonian based from a set t_H_coo-type
    procedure :: init_coo
    procedure :: init_connect    
    procedure :: init_mult_connect_2
end type

abstract interface
    subroutine int_set_from_Hcoo(this,H_coo)
        import t_H_coo,t_H_coo_based,lattice
        class(t_H_coo_based),intent(inout)  :: this
        type(t_H_coo),intent(inout)         :: H_coo
    end subroutine
end interface

contains


subroutine init_coo(this,rowind,colind,val,dim_mode,op_l,op_r,lat,mult_M_single)
    !constructs a Hamiltonian based directly on the coo-arrays, which are moved into the type
    use m_derived_types, only: lattice
    class(t_H_coo_based),intent(inout)        :: this
    real(8),allocatable,intent(inout)   :: val(:)
    integer,allocatable,intent(inout)   :: rowind(:),colind(:)
    integer,intent(in)                  :: dim_mode(2)
    character(len=*),intent(in)         :: op_l         !which order parameters are used at left  side of local Hamiltonian-matrix
    character(len=*),intent(in)         :: op_r         !which order parameters are used at right side of local Hamiltonian-matrix
    type(lattice),intent(in)            :: lat
    integer,intent(in)                  :: mult_M_single !gives the multiple with which the energy_single calculation has to be multiplied (1 for on-site terms, 2 for eg. magnetic exchange)

    type(t_H_coo)           :: H_coo

    Call H_coo%init_coo(rowind,colind,val,dim_mode,op_l,op_r,lat,mult_M_single)
    Call this%init_otherH(H_coo)
    Call this%set_from_Hcoo(H_coo)
end subroutine


subroutine init_connect(this,connect,Hval,Hval_ind,order,lat,mult_M_single)
    use m_derived_types, only: lattice
    class(t_H_coo_based),intent(inout)    :: this
    type(lattice),intent(in)        :: lat
    character(2),intent(in)         :: order
    real(8),intent(in)              :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
    integer,intent(in)              :: Hval_ind(:,:)
    integer,intent(in)              :: connect(:,:)
    integer,intent(in)              :: mult_M_single
    !local
    type(t_H_coo)           :: H_coo

    Call H_coo%init_connect(connect,Hval,Hval_ind,order,lat,mult_M_single)
    Call this%init_otherH(H_coo)
    Call this%set_from_Hcoo(H_coo)
end subroutine 

subroutine init_mult_connect_2(this,connect,Hval,Hval_ind,op_l,op_r,lat,mult_M_single,dim_mode_in)
    !Constructs a Hamiltonian that depends on more than 2 order parameters but only at 2 sites (i.e. some terms are onsite)
    !(example: ME-coupling M_i*E_i*M_j
    use m_derived_types, only: lattice,op_abbrev_to_int
    class(t_H_coo_based),intent(inout)  :: this
    type(lattice),intent(in)            :: lat
    !input Hamiltonian
    real(8),intent(in)                  :: Hval(:)          !values of local Hamiltonian for each line
    integer,intent(in)                  :: Hval_ind(:,:)    !indices in order-parameter space for Hval
    character(len=*),intent(in)         :: op_l             !which order parameters are used at left  side of local Hamiltonian-matrix
    character(len=*),intent(in)         :: op_r             !which order parameters are used at right side of local Hamiltonian-matrix
    integer,intent(in)                  :: connect(:,:)     !lattice sites to be connected (2,Nconnections)
    integer,intent(in)                  :: mult_M_single
    integer,intent(in),optional         :: dim_mode_in(2)   !optional way of putting in dim_mode directly (mainly for custom(not fully unfolded)rankN tensors)
    !local
    type(t_H_coo)           :: H_coo

    Call H_coo%init_mult_connect_2(connect,Hval,Hval_ind,op_l,op_r,lat,mult_M_single,dim_mode_in)
    Call this%init_otherH(H_coo)
    Call this%set_from_Hcoo(H_coo)
end subroutine



end module
