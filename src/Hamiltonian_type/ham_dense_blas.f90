module m_H_dense_blas
#ifdef CPP_BLAS
!Hamiltonian type specifications using dense matrices and no external library
use m_H_deriv, only: t_H
use m_H_dense
use m_derived_types, only: lattice


type,extends(t_H_dense) :: t_H_dense_blas
contains
    procedure :: mult_r,mult_l
end type

private
public t_H,t_H_dense_blas
contains 

subroutine mult_r(this,lat,res,beta)
    !mult
    use m_derived_types, only: lattice
    class(t_H_dense_blas),intent(in)    :: this
    type(lattice), intent(in)           :: lat
    real(8), intent(inout)              :: res(:)   !result matrix-vector product
    real(8),intent(in),optional         :: beta
    ! internal
    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)
    real(8),parameter          :: alpha=1.0d0
    
    Call this%mode_r%get_mode(lat,modes,vec)
    if(size(res)/=this%dimH(1)) STOP "size of vec is wrong"
    if(present(beta))then
        Call dgemm('n','n',this%dimH(1),1,this%dimH(2),alpha,this%H,this%dimH(1),modes,this%dimH(2),beta ,res,this%dimH(1))
    else
        Call dgemm('n','n',this%dimH(1),1,this%dimH(2),alpha,this%H,this%dimH(1),modes,this%dimH(2),0.0d0,res,this%dimH(1))
    endif
    if(allocated(vec)) deallocate(vec)
end subroutine 


subroutine mult_l(this,lat,res,beta)
    use m_derived_types, only: lattice
    class(t_H_dense_blas),intent(in)   :: this
    type(lattice), intent(in)           :: lat
    real(8), intent(inout)              :: res(:)
    real(8),intent(in),optional         :: beta
    ! internal
    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)
    real(8),parameter          :: alpha=1.0d0

    Call this%mode_l%get_mode(lat,modes,vec)
    if(size(res)/=this%dimH(2)) STOP "size of vec is wrong"
    if(present(beta))then
        Call dgemm('t','n',1,this%dimH(2),this%dimH(1),alpha,modes,this%dimh(1),this%H,this%dimH(1),beta ,res,1)
    else
        Call dgemm('t','n',1,this%dimH(2),this%dimH(1),alpha,modes,this%dimh(1),this%H,this%dimH(1),0.0d0,res,1)
    endif
    if(allocated(vec)) deallocate(vec)
end subroutine 

#endif
end module
