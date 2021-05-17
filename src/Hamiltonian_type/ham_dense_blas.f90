module m_H_dense_blas
#ifdef CPP_BLAS
!Hamiltonian type specifications using dense matrices and no external library
use m_H_deriv, only: t_H
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
    
    Call this%mode_r%get_mode(lat,modes,vec)
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

    Call this%mode_l%get_mode(lat,modes,vec)
    if(size(res)/=this%dimH(2)) STOP "size of vec is wrong"
    Call dgemm('t','n',1,this%dimH(2),this%dimH(1),alpha,modes,this%dimh(1),this%H,this%dimH(1),beta,res,1)
    if(allocated(vec)) deallocate(vec)
end subroutine 


subroutine eval_single(this,E,i_m,order,lat)
    use m_derived_types, only: lattice,number_different_order_parameters
    ! input
    class(t_h_dense_blas),intent(in)    :: this
    type(lattice), intent(in)           :: lat
    integer, intent(in)                 :: i_m
    integer, intent(in)                 :: order
!    integer, intent(in)                 :: dim_bnd(2,number_different_order_parameters)    !starting/final index in respective dim_mode of the order parameter (so that energy of single magnetic atom can be be calculated
    ! output
    real(8), intent(out)                :: E
    ! internal
    real(8),pointer                 :: modes_l(:),modes_r(:)
    real(8),allocatable,target      :: vec_l(:),vec_r(:)
    integer                     :: bnd(2)
    integer                     :: size_vec_r
    real(8)                     :: tmp(this%dimH(1))
    real(8),external            :: ddot !blas routine
    real(8),parameter           :: alpha=1.0d0,beta=0.0d0

    ERROR STOP "THIS PROBABLY NO LONGER WORKS WITH THE NEW MODE_L/MODE_R"   !and in general might be much more difficult to implement with eg. rank 4 in M-space only
!    Call lat%point_order(this%op_l,this%dimH(1),modes_l,vec_l)
!    Call lat%point_order_single(this%op_r,i_m,dim_bnd,this%dim_mode(2),modes_r,vec_r,bnd)
!
!    size_vec_r=bnd(2)-bnd(1)+1
!    Call DGEMV('N',this%dimH(1),size_vec_r,alpha,this%H(1,bnd(1)),this%dimH(1),modes_r(1),1,beta,tmp,1)
!    E=ddot(this%dimH(2),modes_l,1,tmp,1)
end subroutine 
#endif
end module
