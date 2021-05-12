module m_H_cusparse
#if defined(CPP_CUDA)
!cuda cusparse implementation
!so far only multiplication to the right works,
!extension to left side can be done by introducing lrev_in/out , lbuffer, etc.
!notice that it requires the transpose matrix, which is supposedly slower, so an overload like 
!  t_H_eigen_mem saving the transpose as well might be sensible for faster evaluation, if that is neccessary

!also single evaluation should be implemented

use m_derived_types, only: lattice,number_different_order_parameters
use m_H_coo_based
use m_cuda_H_interface
use m_cuda_H_vec_interface
use, intrinsic :: iso_c_binding

private
public t_H_cusparse
type(C_PTR),private     ::  handle=c_null_ptr

type,extends(t_H_coo_based) :: t_H_cusparse
    private
    type(C_PTR)     ::  H=c_null_ptr
    type(C_PTR)     ::  rvec_in=c_null_ptr
    type(C_PTR)     ::  rvec_out=c_null_ptr
    type(C_PTR)     ::  buffer=c_null_ptr
contains
!    !necessary t_H routines
    procedure :: eval_single

    procedure :: set_from_Hcoo

    procedure :: add_child 
    procedure :: destroy_child    
    procedure :: mv_child 
    procedure :: copy_child 
    procedure :: bcast_child 

    procedure :: optimize
    procedure :: mult_r,mult_l
    procedure :: mult_l_cont,mult_r_cont
    procedure :: mult_l_disc,mult_r_disc

    !MPI (real implementation might be complicated)
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
        Call cuda_fvec_alloccopy(this%rvec_in,  Hout%rvec_in)
        Call cuda_fvec_alloccopy(this%rvec_out, Hout%rvec_out)
        Call cuda_set_buffer(Hout%buffer,Hout%H,Hout%rvec_in,Hout%rvec_out,handle)
        Call this%copy_deriv(Hout)
    class default
        STOP "Cannot copy t_H_cusparse type to Hamiltonian that is not a class of t_H_cusparse"
    end select
end subroutine

subroutine bcast_child(this,comm)
    use mpi_basic
    use mpi_util,only: bcast
!    use mkl_spblas_util, only: unpack_csr
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
        Call cuda_free_buffer(this%buffer)
        Call cuda_set_buffer(this%buffer,this%H,this%rvec_in,this%rvec_out,handle)
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
        if(c_associated(this%rvec_in)) Call cuda_fvec_destroy(this%rvec_in)
        this%rvec_in=c_null_ptr
        if(c_associated(this%rvec_out)) Call cuda_fvec_destroy(this%rvec_out)
        this%rvec_out=c_null_ptr
        Call cuda_free_buffer(this%buffer)
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
    Call cuda_fvec_init(this%rvec_in,dimH(2))
    Call cuda_fvec_init(this%rvec_out,dimH(2))
    Call cuda_set_buffer(this%buffer,this%H,this%rvec_in,this%rvec_out,handle)
end subroutine 

subroutine eval_single(this,E,i_m,dim_bnd,lat)
    use m_derived_types, only: lattice
    ! input
    class(t_H_cusparse),intent(in)   :: this
    type(lattice), intent(in)       :: lat
    integer, intent(in)             :: i_m
    integer, intent(in)             :: dim_bnd(2,number_different_order_parameters)  !not implemented
    ! output
    real(8), intent(out)            :: E
    ERROR STOP "IMPLEMENT"
end subroutine 

subroutine mult_r(this,lat,res)
    !mult
    use m_derived_types, only: lattice
    class(t_H_cusparse),intent(in)  :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)   !result matrix-vector product
    ! internal
    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)

    Call this%mode_r%get_mode(lat,modes,vec)
    Call cuda_fvec_set(this%rvec_in,modes)
    Call cuda_H_mult_mat_vec(this%H,this%rvec_in, this%rvec_out, this%buffer, handle)
    Call cuda_fvec_get(this%rvec_out,res)
    if(allocated(vec)) deallocate(vec)
end subroutine 

subroutine mult_l(this,lat,res)
    use m_derived_types, only: lattice
    class(t_H_cusparse),intent(in)  :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)
    ERROR STOP "IMPLEMENT"
end subroutine 

subroutine mult_l_cont(this,bnd,vec,res)
    !multiply to the right with a continuous section of the right vector
    class(t_H_cusparse),intent(in)  :: this
    integer,intent(in)              :: bnd(2)
    real(8),intent(in)              :: vec(bnd(2)-bnd(1)+1)
    real(8),intent(inout)           :: res(:)   !result matrix-vector product

    ERROR STOP "IMPLEMENT"
end subroutine 

subroutine mult_r_cont(this,bnd,vec,res)
    !multiply to the right with a continuous section of the right vector
    class(t_H_cusparse),intent(in)  :: this
    integer,intent(in)              :: bnd(2)
    real(8),intent(in)              :: vec(bnd(2)-bnd(1)+1)
    real(8),intent(inout)           :: res(:)   !result matrix-vector product

    ERROR STOP "IMPLEMENT"
end subroutine 

subroutine mult_l_disc(this,N,ind,vec,res)
    !multiply to the right with a discontinuous section of the right vector
    class(t_H_cusparse),intent(in)  :: this
    integer,intent(in)              :: N
    integer,intent(in)              :: ind(N)
    real(8),intent(in)              :: vec(N)
    real(8),intent(inout)           :: res(:)   !result matrix-vector product

    ERROR STOP "IMPLEMENT"
end subroutine 

subroutine mult_r_disc(this,N,ind,vec,res)
    !multiply to the right with a discontinuous section of the right vector
    class(t_H_cusparse),intent(in)  :: this
    integer,intent(in)              :: N
    integer,intent(in)              :: ind(N)
    real(8),intent(in)              :: vec(N)
    real(8),intent(inout)           :: res(:)   !result matrix-vector product

    ERROR STOP "IMPLEMENT"
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
    use mkl_spblas_util, only: unpack_csr
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
    use mkl_spblas_util, only: unpack_csr
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
    use mkl_spblas_util, only: unpack_csr
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
