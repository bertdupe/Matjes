module m_H_dense_blas
#ifdef CPP_BLAS
!Hamiltonian type specifications using dense matrices and no external library
use m_H_combined, only: t_H
use m_H_dense
use m_work_ham_single, only:  work_mode, N_work
use m_derived_types, only: lattice


type,extends(t_H_dense) :: t_H_dense_blas
contains
    procedure :: mult_r,mult_l
end type

private
public t_H,t_H_dense_blas
contains 

subroutine mult_r(this,lat,res,work,alpha,beta)
    !mult
    use m_derived_types, only: lattice
    class(t_H_dense_blas),intent(in)    :: this
    type(lattice), intent(in)           :: lat
    real(8), intent(inout)              :: res(:)   !result matrix-vector product
    type(work_mode),intent(inout)       :: work
    real(8),intent(in),optional         :: alpha
    real(8),intent(in),optional         :: beta
    ! internal
    real(8),pointer ,contiguous         :: modes(:)
    real(8)                    :: alp, bet
    integer                    :: work_size(N_work)

    if(present(alpha))then
        alp=alpha
    else
        alp=1.0d0
    endif
    if(present(beta))then
        bet=beta
    else
        bet=0.0d0
    endif
    Call this%mode_r%get_mode(lat,modes,work,work_size)
    if(size(res)/=this%dimH(1)) STOP "size of vec is wrong"
    Call dgemm('n','n',this%dimH(1),1,this%dimH(2),alp,this%H,this%dimH(1),modes,this%dimH(2),bet,res,this%dimH(1))
    nullify(modes)
    work%offset=work%offset-work_size
end subroutine 


subroutine mult_l(this,lat,res,work,alpha,beta)
    use m_derived_types, only: lattice
    class(t_H_dense_blas),intent(in)    :: this
    type(lattice), intent(in)           :: lat
    real(8), intent(inout)              :: res(:)
    type(work_mode),intent(inout)       :: work
    real(8),intent(in),optional         :: alpha
    real(8),intent(in),optional         :: beta
    ! internal
    real(8),pointer,contiguous          :: modes(:)
    real(8)                    :: alp, bet
    integer                    :: work_size(N_work)

    if(present(alpha))then
        alp=alpha
    else
        alp=1.0d0
    endif
    if(present(beta))then
        bet=beta
    else
        bet=0.0d0
    endif
    Call this%mode_l%get_mode(lat,modes,work,work_size)
    if(size(res)/=this%dimH(2)) STOP "size of vec is wrong"
    Call dgemm('t','n',1,this%dimH(2),this%dimH(1),alp,modes,this%dimh(1),this%H,this%dimH(1),bet ,res,1)
    nullify(modes)
    work%offset=work%offset-work_size
end subroutine 

#endif
end module
