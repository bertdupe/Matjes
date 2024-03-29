module m_H_cusparse
#if defined(CPP_CUDA)
!cuda cusparse implementation
!notice that the mult_l requires the transpose matrix, which is supposedly slower, so an overload like 
!  t_H_eigen_mem saving the transpose as well might be sensible for faster evaluation, if that is neccessary

!single evaluation does not really make sense to implement for cuda, as this is not really good parallelizable and the state memcpy's are probably very prohibitive

use m_derived_types, only: lattice,number_different_order_parameters
use m_H_coo_based
use m_cuda_H_interface
use m_cuda_H_vec_interface
use m_work_ham_single, only:  work_mode, N_work
use, intrinsic :: iso_c_binding

private
public t_H_cusparse
type(C_PTR),private     ::  handle=c_null_ptr

type,extends(t_H_coo_based) :: t_H_cusparse
    private
    type(C_PTR)     ::  H=c_null_ptr
    type(C_PTR)     ::  rvec=c_null_ptr
    type(C_PTR)     ::  lvec=c_null_ptr
    type(C_PTR)     ::  buffer_r=c_null_ptr
    type(C_PTR)     ::  buffer_l=c_null_ptr
contains

    !initialization routine
    procedure :: set_from_Hcoo

    procedure :: optimize
    procedure :: mult_r,mult_l

    !utility routines
    procedure :: add_child 
    procedure :: destroy_child    
    procedure :: mv_child 
    procedure :: copy_child 
    procedure :: bcast_child 

    !MPI NOT IMPLEMENTED, just placeholders(real implementation might be complicated)
    procedure :: send
    procedure :: recv
    procedure :: bcast
    procedure :: distribute 
end type

contains

subroutine mv_child(this,Hout)
    class(t_H_cusparse),intent(inout)   :: this
    class(t_H_base),intent(inout)       :: Hout
    integer         ::  stat
    
    select type(Hout)
    class is(t_H_cusparse)
        Hout%H  =this%H
        Hout%H  =c_null_ptr
    class default
        STOP "Cannot mv t_h_cusparse type with Hamiltonian that is not a class of t_h_cusparse"
    end select
end subroutine

subroutine copy_child(this,Hout)
    class(t_H_cusparse),intent(in)  :: this
    class(t_H_base),intent(inout)   :: Hout
    integer         ::  stat
  
    select type(Hout)
    class is(t_H_cusparse)
        Call cuda_H_copy(this%H,Hout%H)
        !COPY VEC
        Call cuda_fvec_alloccopy(this%rvec,  Hout%rvec)
        Call cuda_fvec_alloccopy(this%lvec, Hout%lvec)
        Call cuda_set_buffer(Hout%buffer_r,Hout%H,logical(.false.,C_BOOL),Hout%rvec,Hout%lvec,handle)
        Call cuda_set_buffer(Hout%buffer_l,Hout%H,logical(.true. ,C_BOOL),Hout%lvec,Hout%rvec,handle)
        Call this%copy_deriv(Hout)
    class default
        STOP "Cannot copy t_H_cusparse type to Hamiltonian that is not a class of t_H_cusparse"
    end select
end subroutine

subroutine bcast_child(this,comm)
    use mpi_basic
    use mpi_util,only: bcast
    class(t_H_cusparse),intent(inout)    ::  this
    type(mpi_type),intent(in)           ::  comm
#ifdef CPP_MPI

    ERROR STOP "IMPLEMENT"
#else
    continue
#endif
end subroutine 

subroutine add_child(this,H_in)
    class(t_H_cusparse),intent(inout)   :: this
    class(t_H_base),intent(in)          :: H_in
    
    type(C_PTR)       :: tmp_H

    tmp_H=c_null_ptr
    select type(H_in)
    class is(t_H_cusparse)
        Call cuda_H_add(this%H,H_in%H,tmp_H,handle)
        Call cuda_H_destroy(this%H)
        this%H=tmp_H
        Call cuda_free_buffer(this%buffer_r)
        Call cuda_free_buffer(this%buffer_l)
        Call cuda_set_buffer(this%buffer_r,this%H,logical(.false.,C_BOOL),this%rvec,this%lvec,handle)
        Call cuda_set_buffer(this%buffer_l,this%H,logical(.true. ,C_BOOL),this%lvec,this%rvec,handle)
    class default
        STOP "Cannot add t_h_cusparse type with Hamiltonian that is not a class of t_h_cusparse"
    end select
end subroutine 

subroutine destroy_child(this)
    class(t_H_cusparse),intent(inout)   :: this
    integer     ::  stat

    if(this%is_set())then
        if(c_associated(this%H)) Call cuda_H_destroy(this%H)
        this%H=c_null_ptr
        if(c_associated(this%rvec)) Call cuda_fvec_destroy(this%rvec)
        this%rvec=c_null_ptr
        if(c_associated(this%lvec)) Call cuda_fvec_destroy(this%lvec)
        this%lvec=c_null_ptr
        Call cuda_free_buffer(this%buffer_r)
        Call cuda_free_buffer(this%buffer_l)
    endif
end subroutine

subroutine set_from_Hcoo(this,H_coo)
    class(t_H_cusparse),intent(inout)   :: this
    type(t_H_coo),intent(inout)         :: H_coo
    !local
    integer                 :: nnz
    integer(C_int)          :: stat
    integer(C_int)          :: dimH(2)
    real(C_DOUBLE),allocatable  :: val(:)
    integer,allocatable         :: rowind(:),colind(:)

    Call this%init_otherH(H_coo)
    Call H_coo%pop_par(this%dimH,nnz,val,rowind,colind)
    dimH=this%dimH
    if(.not.c_associated(handle)) Call cuda_create_handle(handle)
    Call cuda_H_init(nnz,dimH,rowind,colind,val,this%H,handle)
    Call cuda_fvec_init(this%lvec,dimH(1))
    Call cuda_fvec_init(this%rvec,dimH(2))
    Call cuda_set_buffer(this%buffer_r,this%H,logical(.false.,C_BOOL),this%rvec,this%lvec,handle)
    Call cuda_set_buffer(this%buffer_l,this%H,logical(.true. ,C_BOOL),this%lvec,this%rvec,handle)
end subroutine 

!subroutine eval_single(this,E,i_m,order,lat,work)
!    use m_work_ham_single
!    use m_derived_types, only: lattice
!    use,intrinsic :: ISO_FORTRAN_ENV, only: error_unit
!    ! input
!    class(t_H_cusparse),intent(in)     :: this
!    type(lattice), intent(in)       :: lat
!    integer, intent(in)             :: i_m
!    integer,intent(in)              :: order
!    ! output
!    real(8), intent(out)            :: E
!    type(work_ham_single),intent(inout) :: work
!
!    write(error_unit,"(///A)") "The eval_single routine is not implemented for cuda Hamiltonians (t_H_cusparse)"
!    write(error_unit,"(A)") "There is no reasonable way to parallelize this evaluation, for faster Monte-Carlo sampling the entire MCstep routine should be parallelized in CUDA"
!    write(error_unit,"(A)") "USE A DIFFERENT Hamiltonian_mode"
!    STOP
!end subroutine 

subroutine mult_r(this,lat,res,work,alpha,beta)
    use m_derived_types, only: lattice
    class(t_H_cusparse),intent(in)  :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)   !result matrix-vector product
    type(work_mode),intent(inout)   :: work
    real(8),intent(in),optional     :: alpha
    real(8),intent(in),optional     :: beta
    ! internal
    real(8),pointer,contiguous :: modes(:)
    real(8)                    :: alp, bet
    integer                    :: work_size(N_work)

    if(present(alpha))then
        alp=alpha
    else
        alp=1.0d0
    endif
    if(present(beta))then
        bet=beta
        Call cuda_fvec_set(this%lvec,res)
    else
        bet=0.0d0
    endif
    Call this%mode_r%get_mode(lat,modes,work,work_size)
    Call cuda_fvec_set(this%rvec,modes)
    Call cuda_H_mult_r(this%H,this%rvec, this%lvec, alp, bet, this%buffer_r, handle)
    Call cuda_fvec_get(this%lvec,res)
    nullify(modes)
    work%offset=work%offset-work_size
end subroutine 

subroutine mult_l(this,lat,res,work,alpha,beta)
    use m_derived_types, only: lattice
    class(t_H_cusparse),intent(in)  :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)
    type(work_mode),intent(inout)   :: work
    real(8),intent(in),optional     :: alpha
    real(8),intent(in),optional     :: beta
    ! internal
    real(8),pointer,contiguous :: modes(:)
    real(8)                    :: alp, bet
    integer                    :: work_size(N_work)

    if(present(alpha))then
        alp=alpha
    else
        alp=1.0d0
    endif
    if(present(beta))then
        bet=beta
        Call cuda_fvec_set(this%rvec,res)
    else
        bet=0.0d0
    endif
    Call this%mode_l%get_mode(lat,modes,work,work_size)
    Call cuda_fvec_set(this%lvec,modes)
    Call cuda_H_mult_l(this%H,this%lvec, this%rvec, alp, bet, this%buffer_l, handle)
    Call cuda_fvec_get(this%rvec,res)
    nullify(modes)
    work%offset=work%offset-work_size
end subroutine 

subroutine optimize(this)
    class(t_H_cusparse),intent(inout)   :: this
   
    
    !write(*,'(///A///)') "MAYBE THERE IS SOMETHING TO DO TO OPTIMIZE THE CUSPARSE MATRIX MULTIPLICATION"
    !ERROR STOP "IMPLEMENT"
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!            MPI ROUTINES           !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine send(this,ithread,tag,com)
    use mpi_basic                
    class(t_H_cusparse),intent(in)   :: this
    integer,intent(in)              :: ithread
    integer,intent(in)              :: tag
    integer,intent(in)              :: com

#ifdef CPP_MPI
    Call this%send_base(ithread,tag,com)

    ERROR STOP "IMPLEMENT MPI for t_H_cusparse, if really necessary"
#else
    continue
#endif
end subroutine

subroutine recv(this,ithread,tag,com)
    use mpi_basic                
    class(t_H_cusparse),intent(inout)   :: this
    integer,intent(in)              :: ithread
    integer,intent(in)              :: tag
    integer,intent(in)              :: com
#ifdef CPP_MPI

    Call this%recv_base(ithread,tag,com)
    ERROR STOP "IMPLEMENT MPI for t_H_cusparse, if really necessary"
#else
    continue
#endif
end subroutine

subroutine bcast(this,comm)
    use mpi_basic
    use mpi_util,only: bcast_util => bcast
    class(t_H_cusparse),intent(inout)    ::  this
    type(mpi_type),intent(in)           ::  comm
#ifdef CPP_MPI
    Call this%bcast_base(comm)
    ERROR STOP "IMPLEMENT MPI for t_H_cusparse, if really necessary"
#else
    continue
#endif
end subroutine 

subroutine distribute(this,comm)
    use mpi_basic                
    use mpi_util!,only: bcast_util => bcast
    class(t_H_cusparse),intent(inout)        ::  this
    type(mpi_type),intent(in)       ::  comm
#ifdef CPP_MPI
    Call this%bcast_base(comm)
    ERROR STOP "IMPLEMENT MPI for t_H_cusparse, if really necessary"
#else
    continue
#endif
end subroutine 

#endif
end module
