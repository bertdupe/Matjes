module m_H_eigen
#ifdef CPP_EIGEN_H
!Hamiltonian type specifications using dense matrices and no external library
use m_H_type
use m_H_coo
use m_derived_types, only: lattice
use m_eigen_H_interface

type,extends(t_H) :: t_H_eigen
    type(C_PTR)     ::  H=c_null_ptr
!some pointer to Hamiltonian
contains
    !necessary t_H routines
    procedure :: eval_single
    procedure :: init      
    procedure :: init_1    
    procedure :: init_mult_2   

    procedure :: add_child 
    procedure :: destroy_child    
    procedure :: copy_child 

    procedure :: optimize
    procedure :: mult_r,mult_l
end type

interface t_H_mkl_csr
    procedure :: dummy_constructor
end interface 

private
public t_H,t_H_eigen
contains 

type(t_H_eigen) function dummy_constructor()
    !might want some initialization for H and descr, but should work without
    !continue 
end function 

subroutine mult_r(this,lat,res)
    !mult
    use m_derived_types, only: lattice
    class(t_H_eigen),intent(in)     :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)   !result matrix-vector product
    ! internal
    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)

    Call lat%point_order(this%op_r,this%dimH(2),modes,vec)
    Call eigen_H_mult_mat_vec(this%H,modes,res)
end subroutine 


subroutine mult_l(this,lat,res)
    use m_derived_types, only: lattice
    class(t_H_eigen),intent(in)     :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)
    ! internal
    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)

    Call lat%point_order(this%op_l,this%dimH(1),modes,vec)
    Call eigen_H_mult_vec_mat(this%H,modes,res)
end subroutine 


subroutine optimize(this)
    class(t_H_eigen),intent(inout)   :: this

    !THERE MIGHT BE SOMETHING TO DO HERE
    continue 
end subroutine

subroutine init_1(this,line,Hval,Hval_ind,order,lat)
    use m_derived_types, only: lattice
    class(t_H_eigen),intent(inout)  :: this

    type(lattice),intent(in)        :: lat
    integer,intent(in)              :: order(2)
    real(8),intent(in)              :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
    integer,intent(in)              :: Hval_ind(:,:)
    integer,intent(in)              :: line(:,:)
    !local
    type(t_H_coo)           :: H_coo

    if(this%is_set()) STOP "cannot set hamiltonian as it is already set"
    Call H_coo%init_1(line,Hval,Hval_ind,order,lat)
    Call set_from_Hcoo(this,H_coo,lat)
end subroutine 

subroutine init_mult_2(this,connect,Hval,Hval_ind,op_l,op_r,lat)
    use m_derived_types, only: lattice
    class(t_H_eigen),intent(inout)  :: this

    type(lattice),intent(in)        :: lat
    integer,intent(in)              :: op_l(:),op_r(:)
    real(8),intent(in)              :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
    integer,intent(in)              :: Hval_ind(:,:)
    integer,intent(in)              :: connect(:,:)
    !local
    type(t_H_coo)           :: H_coo

    if(this%is_set()) STOP "cannot set hamiltonian as it is already set"
    Call H_coo%init_mult_2(connect,Hval,Hval_ind,op_l,op_r,lat)
    Call set_from_Hcoo(this,H_coo,lat)
end subroutine 

subroutine copy_child(this,Hout)
    class(t_H_eigen),intent(in)     :: this
    class(t_H),intent(inout)        :: Hout
    
    select type(Hout)
    class is(t_H_eigen)
        Call eigen_H_copy(this%H,Hout%H) 
    class default
        STOP "Cannot copy t_H_eigen type with Hamiltonian that is not a class of t_H_eigen"
    end select
end subroutine

subroutine add_child(this,H_in)
    class(t_H_eigen),intent(inout)    :: this
    class(t_H),intent(in)             :: H_in

    select type(H_in)
    class is(t_H_eigen)
        Call eigen_H_add(this%H,H_in%H)
    class default
        STOP "Cannot add t_H_eigen type with Hamiltonian that is not a class of t_H_eigen"
    end select
end subroutine 

subroutine destroy_child(this)
    class(t_H_eigen),intent(inout)    :: this

    if(this%is_set())then
        Call eigen_H_destroy(this%H)
    endif
end subroutine

subroutine init(this,energy_in,lat)
    use m_derived_types, only: operator_real_order_N,lattice
    class(t_H_eigen),intent(inout)      :: this
    type(operator_real_order_N)         :: energy_in
    type(lattice),intent(in)            :: lat
    !local
    type(t_H_coo)           :: H_coo

    if(this%is_set()) STOP "cannot set hamiltonian as it is already set"
    Call H_coo%init(energy_in,lat)
    Call set_from_Hcoo(this,H_coo,lat)
end subroutine 

subroutine set_from_Hcoo(this,H_coo,lat)
    type(t_H_eigen),intent(inout)       :: this
    type(t_H_coo),intent(inout)         :: H_coo
    type(lattice),intent(in)            :: lat
    !local
    integer                 :: nnz,i
    real(8),allocatable     :: val(:)
    integer,allocatable     :: rowind(:),colind(:)
    type(C_PTR)     ::  H

    Call H_coo%pop_par(this%dimH,nnz,val,rowind,colind)
    colind=colind-1;rowind=rowind-1
    Call eigen_H_init(nnz,this%dimH,rowind,colind,val,this%H)
    Call this%init_base(lat,H_coo%op_l,H_coo%op_r)
end subroutine 

subroutine eval_single(this,E,i_m,lat)
    use m_derived_types, only: lattice
    ! input
    class(t_H_eigen),intent(in)     :: this
    type(lattice), intent(in)       :: lat
    integer, intent(in)             :: i_m
    ! output
    real(kind=8), intent(out)       :: E
    ! internal
    real(8),pointer                 :: modes_l(:),modes_r(:)
    real(8),allocatable,target      :: vec_l(:),vec_r(:)
    real(8)                         :: tmp(this%dimH(2))

    Call lat%point_order(this%op_l,this%dimH(1),modes_l,vec_l)
    Call lat%point_order(this%op_r,this%dimH(2),modes_r,vec_r)

    ERROR STOP "NOT IMPLEMENTED" 
end subroutine 
#endif 
end module
