module m_H_eigen
#ifdef CPP_EIGEN
!Hamiltonian type specifications using Eigen to save and evaluate the Hamiltonian
use m_derived_types, only: lattice, number_different_order_parameters
use m_eigen_H_interface
use m_H_coo_based
use m_mode_construction_rankN_sparse_col
use m_work_ham_single
use,intrinsic :: ISO_FORTRAN_ENV, only: error_unit

type,extends(t_H_coo_based) :: t_H_eigen
    type(C_PTR)     ::  H=c_null_ptr  !pointer to sparse matrix c-object 

    !pointers to Hamiltonian data handled by eigen (col major format) (zero based)
    real(C_DOUBLE),pointer  :: H_val(:)
    integer(C_INT),pointer  :: H_outer(:),H_inner(:)
    integer                 :: nnz=0
contains
    !necessary t_H routines
    procedure :: set_from_Hcoo

    procedure :: add_child 
    procedure :: destroy_child
    procedure :: copy_child 
    procedure :: optimize

    procedure :: mult_l, mult_r

    procedure :: mult_l_disc, mult_r_disc

    procedure :: set_work_single
    procedure :: get_work_size_single

    procedure :: set_auxiliaries

    !MPI
    procedure :: send
    procedure :: recv
    procedure :: distribute
    procedure :: bcast
end type

interface t_H_mkl_csr
    procedure :: dummy_constructor
end interface 

private
public t_H,t_H_eigen
contains 

subroutine set_work_single(this,work,order)
    class(t_H_eigen),intent(inout)          :: this
    class(work_ham_single),intent(inout)    :: work 
    integer,intent(in)                      :: order
    integer     :: sizes(2)
    integer     :: dim_mode

    if(.not.this%is_set()) ERROR STOP "cannot set work size of hamiltonian if it is not set"
    if(this%col_max==0) ERROR STOP "cannot set work size of t_H_mkl_csr if col_max==0"
    Call this%get_work_size_single(sizes)
    Call work%set(sizes)
end subroutine

subroutine get_work_size_single(this,sizes)
    class(t_H_eigen),intent(in) :: this
    integer,intent(out)         :: sizes(2)

    Call work_size_single(maxval(this%dim_r_single),this%col_max,sizes)
end subroutine

subroutine mult_l_disc(this,i_m,lat,N,ind_out,vec,ind_sum,ind_Mult,mat_mult,vec_mult)
    !Calculates the entries of the left vector * matrix product for the indices ind_out of the result vector
    USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
    ! input
    class(t_H_eigen), intent(in)        :: this
    type(lattice), intent(in)           :: lat
    integer, intent(in)                 :: i_m          !index of the comp's right mode in the inner dim_mode
    integer, intent(in)                 :: N            !number of indices to calculated
    integer, intent(in)                 :: ind_out(N)   !indices to be calculated
    ! output
    real(8),intent(out)                 :: vec(N)       !vector.matrix result
    !temporary data
    integer,intent(inout)               :: ind_sum (N+1)                !ind_mult index where the a vec-entry start and end
    integer,intent(inout)               :: ind_mult(N*this%col_max)     !indices of the left array which have non-vanishing contributions to get the ind entries of the vec/mat product
    real(8),intent(inout)               :: mat_mult(N*this%col_max)     !matrix entries corresponding to ind_mult
    real(8),intent(inout)               :: vec_mult(N*this%col_max)     !values of discontiguous left mode which has to be evaluated (indices of ind_mult)

    !some local indices/ loop variables
    integer ::  i,j, i_outer,ii

    !get matrix indices and values whose vec.mat product constitute the output vector
    ind_sum(1)=0
    ii=0
    do i=1,N
        i_outer=ind_out(i)
        do j=this%H_outer(i_outer)+1,this%H_outer(i_outer+1)
            ii=ii+1
            mat_mult(ii)=this%H_val(j)
            ind_mult(ii)=this%H_inner(j)
        enddo
        ind_sum(i+1)=ii
    enddo
    ind_mult=ind_mult+1    !zero based index from c++

    !get the left vector values which are multiplied with the matrix
    Call this%mode_l%get_mode_disc(lat,ii,ind_mult(:ii),vec_mult(:ii))

    !multipy the matrix and left vector entries and sum together to the respective output entry
    vec_mult(:ii)=vec_mult(:ii)*mat_mult(:ii)
    do i=1,N
        vec(i)=sum(vec_mult(ind_sum(i)+1:ind_sum(i+1)))
    enddo
end subroutine 

subroutine mult_r_disc(this,i_m,lat,N,ind_out,vec,ind_sum,ind_Mult,mat_mult,vec_mult)
    !Calculates the entries of the matrix * right vector product for the indices ind_out of the result vector
    ! input
    class(t_H_eigen), intent(in)        :: this
    type(lattice), intent(in)           :: lat
    integer, intent(in)                 :: i_m
    integer, intent(in)                 :: N
    integer, intent(in)                 :: ind_out(N)
    ! output
    real(8),intent(out)                 :: vec(N)
    !temporary data
    integer,intent(inout)               :: ind_sum (N+1)
    integer,intent(inout)               :: ind_mult(N*this%row_max)
    real(8),intent(inout)               :: mat_mult(N*this%row_max)
    real(8),intent(inout)               :: vec_mult(N*this%row_max)

    write(error_unit,'(//A)') "Trying to call mult_l_dist from type t_H_eigen."
    write(error_unit,'(A)')   "This can not be done efficiently without the transpose as implemented in t_H_eigen_mem."
    write(error_unit,'(A)')   "Please choose the Hamiltonian implementation including the transpose (Hamiltonian_mode)."
    ERROR STOP
end subroutine 

type(t_H_eigen) function dummy_constructor()
    !might want some initialization for H and descr, but should work without

    !continue 
end function 

subroutine mult_r(this,lat,res,alpha,beta)
    class(t_H_eigen),intent(in)     :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)   !result matrix-vector product
    real(8),intent(in),optional     :: alpha
    real(8),intent(in),optional     :: beta
    ! internal
    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)
    real(8)                    :: alp, bet

    if(present(alpha))then
        alp=alpha
    else
        alp=1.0d0
    endif
    if(present(beta))then
        bet=beta
    else
        res=0.0d0  !prevents uninitialized runtime warnings
        bet=0.0d0
    endif
    Call this%mode_r%get_mode(lat,modes,vec)
    Call eigen_H_mult_r(this%H,modes,res,alp,bet)
    if(allocated(vec)) deallocate(vec)
end subroutine 

subroutine mult_r_ind(this,lat,N,ind_out,vec_out)
    class(t_H_eigen),intent(in)     :: this
    type(lattice),intent(in)        :: lat
    integer,intent(in)              :: N
    integer,intent(in)              :: ind_out(N)
    real(8),intent(out)             :: vec_out(N)

    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)

    Call this%mode_r%get_mode(lat,modes,vec)
    Call eigen_H_mult_r_ind(this%H,modes,N,ind_out,vec_out)
end subroutine

subroutine mult_l_ind(this,lat,N,ind_out,vec_out)
    class(t_H_eigen),intent(in)     :: this
    type(lattice),intent(in)        :: lat
    integer,intent(in)              :: N
    integer,intent(in)              :: ind_out(N)
    real(8),intent(out)             :: vec_out(N)

    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)

    Call this%mode_l%get_mode(lat,modes,vec)
    Call eigen_H_mult_l_ind(this%H,modes,N,ind_out,vec_out)
end subroutine

subroutine mult_l(this,lat,res,alpha,beta)
    use m_derived_types, only: lattice
    class(t_H_eigen),intent(in)     :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)
    real(8),intent(in),optional     :: alpha
    real(8),intent(in),optional     :: beta
    ! internal
    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)
    real(8)                    :: alp, bet

    if(present(alpha))then
        alp=alpha
    else
        alp=1.0d0
    endif
    if(present(beta))then
        bet=beta
    else
        res=0.0d0  !prevents uninitialized runtime warnings
        bet=0.0d0
    endif
    Call this%mode_l%get_mode(lat,modes,vec)
    Call eigen_H_mult_l(this%H,modes,res,alp,bet)
    if(allocated(vec)) deallocate(vec)
end subroutine 

subroutine optimize(this)
    class(t_H_eigen),intent(inout)   :: this

    !THERE MIGHT BE SOMETHING TO DO HERE
    continue 
end subroutine

subroutine copy_child(this,Hout)
    class(t_H_eigen),intent(in)     :: this
    class(t_H_base),intent(inout)   :: Hout
    
    select type(Hout)
    class is(t_H_eigen)
        Call eigen_H_copy(this%H,Hout%H) 
        Call set_auxiliaries(Hout)
    class default
        STOP "Cannot copy t_H_eigen type with Hamiltonian that is not a class of t_H_eigen"
    end select
end subroutine

subroutine add_child(this,H_in)
    class(t_H_eigen),intent(inout)    :: this
    class(t_H_base),intent(in)        :: H_in

    select type(H_in)
    class is(t_H_eigen)
        Call eigen_H_add(this%H,H_in%H)
        Call this%set_auxiliaries()
    class default
        STOP "Cannot add t_H_eigen type with Hamiltonian that is not a class of t_H_eigen"
    end select
end subroutine 

subroutine destroy_child(this)
    class(t_H_eigen),intent(inout)    :: this

    if(this%is_set())then
        Call eigen_H_destroy(this%H)
        nullify(this%H_val,this%H_outer,this%H_inner)
        this%col_max=0
        this%nnz=0
    endif
end subroutine

subroutine set_from_Hcoo(this,H_coo)
    class(t_H_eigen),intent(inout)  :: this
    type(t_H_coo),intent(inout)     :: H_coo
    !local
    integer                 :: nnz
    real(8),allocatable     :: val(:)
    integer,allocatable     :: rowind(:),colind(:)

    Call H_coo%pop_par(this%dimH,nnz,val,rowind,colind)
    colind=colind-1;rowind=rowind-1
    Call eigen_H_init(nnz,this%dimH,rowind,colind,val,this%H)
    Call this%set_auxiliaries()
end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!            MPI ROUTINES           !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine send(this,ithread,tag,com)
    use mpi_basic                
    class(t_H_eigen),intent(in)     :: this
    integer,intent(in)              :: ithread
    integer,intent(in)              :: tag
    integer,intent(in)              :: com

#ifdef CPP_MPI
    Call this%send_base(ithread,tag,com)
    Call eigen_H_send(ithread,tag,this%H,com) 
#else
    continue
#endif
end subroutine

subroutine recv(this,ithread,tag,com)
    use mpi_basic                
    class(t_H_eigen),intent(inout)   :: this
    integer,intent(in)              :: ithread
    integer,intent(in)              :: tag
    integer,intent(in)              :: com

#ifdef CPP_MPI
    Call this%recv_base(ithread,tag,com)
    Call eigen_H_recv(ithread,tag,this%H,com) 
    Call this%set_auxiliaries()
#else
    continue
#endif
end subroutine

subroutine bcast(this,comm)
    use mpi_basic                
    class(t_H_eigen),intent(inout)  ::  this
    type(mpi_type),intent(in)       ::  comm
#ifdef CPP_MPI
    Call this%bcast_base(comm)
    Call eigen_H_bcast(comm%id,comm%mas,comm%ismas,this%H,comm%com) 
    if(.not.comm%ismas) Call this%set_auxiliaries()
#else
    continue
#endif
end subroutine 

subroutine distribute(this,comm)
    use mpi_basic                
    class(t_H_eigen),intent(inout)        ::  this
    type(mpi_type),intent(in)       ::  comm

#ifdef CPP_MPI
    Call this%bcast_base(comm)
    Call eigen_H_distribute(comm%id,comm%mas,comm%ismas,this%H,comm%com)
    if(.not.comm%ismas) Call this%set_auxiliaries()
#else
    continue
#endif
end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           HELPER FUNCTIONS                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_auxiliaries(this)
    class(t_H_eigen),intent(inout)    :: this

    Call this%set_deriv()
    Call set_H_ptr(this)
    Call set_col_max(this)
end subroutine

subroutine set_H_ptr(this)
    class(t_H_eigen),intent(inout)    :: this

    type(C_PTR)     :: col,row,val
    integer         :: dimH(2)

    nullify(this%H_inner,this%H_outer,this%H_val)
    Call eigen_get_dat(this%H,this%nnz,dimH,col,row,val)
    CAll C_F_POINTER(col, this%H_outer, [dimH(2)+1])
    CAll C_F_POINTER(row, this%H_inner, [this%nnz]) 
    CAll C_F_POINTER(val, this%H_val,   [this%nnz]) 
end subroutine

subroutine set_col_max(this)
    class(t_H_eigen),intent(inout)    :: this
    !variable to get multiplication for single evaluation (estimate sizes...)
    integer,allocatable               :: tmp(:)
    integer                           :: i

    allocate(tmp(size(this%H_outer)-1))
    do i=1,size(tmp)
        tmp(i)=this%H_outer(i+1)-this%H_outer(i)
    enddo
    this%col_max=maxval(tmp)
end subroutine
#endif 
end module
