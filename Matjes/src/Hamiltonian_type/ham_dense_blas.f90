module m_H_dense_blas
#ifdef CPP_BLAS
!Hamiltonian type specifications using dense matrices and no external library
use m_H_type
use m_H_dense
use m_derived_types, only: lattice


type,extends(t_H_dense) :: t_H_dense_blas
contains
    procedure :: eval_single
    procedure :: mult_r,mult_l
end type

private
public t_H,t_H_dense_blas
contains 

subroutine mult_r(this,lat,res)
    !mult
    use m_derived_types, only: lattice
    class(t_H_dense_blas),intent(in)   :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)   !result matrix-vector product
    ! internal
    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)
    real(8),parameter          :: alpha=1.0d0,beta=0.0d0
    

    Call lat%point_order(this%op_r,this%dimH(2),modes,vec)
    if(size(res)/=this%dimH(1)) STOP "size of vec is wrong"

    Call dgemm('n','n',this%dimH(1),1,this%dimH(2),alpha,this%H,this%dimH(1),modes,this%dimH(2),beta,res,this%dimH(1))

    if(allocated(vec)) deallocate(vec)
end subroutine 


subroutine mult_l(this,lat,res)
    use m_derived_types, only: lattice
    class(t_H_dense_blas),intent(in)   :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)
    ! internal
    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)
    real(8),parameter          :: alpha=1.0d0,beta=0.0d0

    Call lat%point_order(this%op_l,this%dimH(1),modes,vec)

    if(size(res)/=this%dimH(2)) STOP "size of vec is wrong"

    Call dgemm('t','n',1,this%dimH(2),this%dimH(1),alpha,modes,this%dimh(1),this%H,this%dimH(1),beta,res,1)
    if(allocated(vec)) deallocate(vec)
    
end subroutine 


subroutine eval_single(this,E,i_m,lat)
    use m_derived_types, only: lattice
    ! input
    class(t_h_dense_blas),intent(in)    :: this
    type(lattice), intent(in)           :: lat
    integer, intent(in)                 :: i_m
    ! output
    real(8), intent(out)                :: E
    ! internal
    real(8),pointer                 :: modes_l(:),modes_r(:)
    real(8),allocatable,target      :: vec_l(:),vec_r(:)
    real(8)                     :: tmp(this%dimH(1))
    integer                     :: ind
    real(8),external            :: ddot !blas routine
    real(8),parameter           :: alpha=1.0d0,beta=0.0d0

    Call lat%point_order(this%op_l,this%dimH(1),modes_l,vec_l)
    Call lat%point_order(this%op_r,this%dimH(2),modes_r,vec_r)

    ind=1+(i_m-1)*this%dim_mode(2)
    Call DGEMV('N',this%dimH(1),this%dim_mode(2),alpha,this%H(1,ind),this%dimH(1),modes_r(ind),1,beta,tmp,1)
    E=ddot(this%dimH(2),modes_l,1,tmp,1)
end subroutine 
#endif
end module
