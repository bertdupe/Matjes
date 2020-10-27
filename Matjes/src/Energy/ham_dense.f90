module m_H_dense
!Hamiltonian type specifications using dense matrices and no external library
use m_H_type
use m_H_coo
use m_derived_types, only: lattice


type,extends(t_H) :: t_H_dense
    private
    real(8),allocatable   :: H(:,:)
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

private
public t_H,t_H_dense
contains 

subroutine mult_r(this,lat,res)
    !mult
    use m_derived_types, only: lattice
    class(t_h_dense),intent(in)   :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)   !result matrix-vector product
    ! internal
    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)

    Call lat%point_order(this%op_r,this%dimH(2),modes,vec)

    if(size(res)/=this%dimH(1)) STOP "size of vec is wrong"
    res=matmul(this%H,modes)
    if(allocated(vec)) deallocate(vec)
end subroutine 


subroutine mult_l(this,lat,res)
    use m_derived_types, only: lattice
    class(t_h_dense),intent(in)   :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)
    ! internal
    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)

    Call lat%point_order(this%op_l,this%dimH(1),modes,vec)

    if(size(res)/=this%dimH(2)) STOP "size of vec is wrong"
    res=matmul(modes,this%H)
    if(allocated(vec)) deallocate(vec)
    
end subroutine 


subroutine optimize(this)
    class(t_h_dense),intent(inout)   :: this

    !nothing to optimize here
    continue 
end subroutine

subroutine init_1(this,line,Hval,Hval_ind,order,lat)
    use m_derived_types, only: lattice
    class(t_h_dense),intent(inout)    :: this

    type(lattice),intent(in)        :: lat
    integer,intent(in)              :: order(2)
    real(8),intent(in)              :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
    integer,intent(in)              :: Hval_ind(:,:)
    integer,intent(in)              :: line(:,:)

    !local
    type(t_H_coo)           :: H_coo

    if(this%is_set()) STOP "cannot set hamiltonian as it is already set"
    Call H_coo%init_1(line,Hval,Hval_ind,order,lat)
    Call set_from_Hcoo(this,H_coo)
end subroutine 


subroutine init_mult_2(this,connect,Hval,Hval_ind,op_l,op_r,lat)
    use m_derived_types, only: lattice
    class(t_h_dense),intent(inout)    :: this

    type(lattice),intent(in)        :: lat
    integer,intent(in)              :: op_l(:),op_r(:)
    real(8),intent(in)              :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
    integer,intent(in)              :: Hval_ind(:,:)
    integer,intent(in)              :: connect(:,:)

    !local
    type(t_H_coo)           :: H_coo

    if(this%is_set()) STOP "cannot set hamiltonian as it is already set"
    Call H_coo%init_mult_2(connect,Hval,Hval_ind,op_l,op_r,lat)
    Call set_from_Hcoo(this,H_coo)
end subroutine 


subroutine copy_child(this,Hout)
    class(t_h_dense),intent(in)   :: this
    class(t_H),intent(inout)        :: Hout
    
    select type(Hout)
    class is(t_h_dense)
        allocate(Hout%H,source=this%H)
    class default
        STOP "Cannot copy t_h_dense type with Hamiltonian that is not a class of t_h_dense"
    end select
end subroutine

subroutine add_child(this,H_in)
    class(t_h_dense),intent(inout)    :: this
    class(t_H),intent(in)             :: H_in

    select type(H_in)
    class is(t_h_dense)
        this%H=this%H+H_in%H
    class default
        STOP "Cannot add t_h_dense type with Hamiltonian that is not a class of t_h_dense"
    end select

end subroutine 

subroutine destroy_child(this)
    class(t_h_dense),intent(inout)    :: this

    if(this%is_set())then
        deallocate(this%H)
    endif
end subroutine

subroutine init(this,energy_in,lat)
    use m_derived_types, only: operator_real_order_N,lattice
    class(t_h_dense),intent(inout)    :: this
    type(operator_real_order_N)         :: energy_in
    type(lattice),intent(in)            :: lat

    !local
    type(t_H_coo)           :: H_coo

    if(this%is_set()) STOP "cannot set hamiltonian as it is already set"
    Call H_coo%init(energy_in,lat)
    Call set_from_Hcoo(this,H_coo)

end subroutine 

subroutine set_from_Hcoo(this,H_coo)
    type(t_H_coo),intent(inout)         :: H_coo
    type(t_h_dense),intent(inout)     :: this

    !local
    integer                 :: nnz,i
    real(8),allocatable     :: val(:)
    integer,allocatable     :: rowind(:),colind(:)

    Call H_coo%pop_par(this%dimH,nnz,val,rowind,colind)
    allocate(this%H(this%dimH(1),this%dimH(2)),source=0.0d0)
    do i=1,nnz
        this%H(rowind(i),colind(i))=val(i)
    enddo

    allocate(this%op_l,source=H_coo%op_l)
    allocate(this%op_r,source=H_coo%op_r)
    Call this%set_prepared(.true.)
end subroutine 


subroutine eval_single(this,E,i_m,lat)
    use m_derived_types, only: lattice
    ! input
    class(t_h_dense),intent(in)    :: this
    type(lattice), intent(in)       :: lat
    integer, intent(in)             :: i_m
    ! output
    real(kind=8), intent(out)       :: E
    ! internal
    real(8)                         :: tmp(this%dimH(1))
    real(8),pointer             :: modes_l(:)
    real(8),allocatable,target  :: vec_l(:)
    integer             :: dim_modes(2)

    Call lat%point_order(this%op_l,this%dimH(1),modes_l,vec_l)
    Call this%mult_r(lat,tmp)

    STOP "UPDATE EVAL_SINGLE FOR HIGHER RANKS, AND ALSO DIM_MODES SEEMS BROKEN...ALSO NOT ADJUSTED FOR t_h_dense at all" !there is some function to get those already
    ! try sparse matrix product to substitute sparse matrix times sparse vector product
    E=dot_product(modes_l((i_m-1)*dim_modes(1)+1:i_m*dim_modes(1)),tmp((i_m-1)*dim_modes(1)+1:i_m*dim_modes(1)))
end subroutine 

end module
