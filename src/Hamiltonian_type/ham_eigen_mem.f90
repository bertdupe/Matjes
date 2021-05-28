module m_H_eigen_mem
#ifdef CPP_EIGEN
!Hamiltonian type specifications extending the Eigen version to also save the transpose 
!in order to make some operations faster
use m_derived_types, only: lattice, number_different_order_parameters
use m_eigen_H_interface
use m_H_coo_based
use m_H_eigen, only: t_h_eigen
use m_mode_construction_rankN_sparse_col
use m_work_ham_single

type,extends(t_H_eigen) :: t_H_eigen_mem
    type(C_PTR)             ::  HT=c_null_ptr  !transpose Hamiltonian in eigen

    !pointers to Hamiltonian data handled by eigen (col major format) (zero based)
    real(C_DOUBLE),pointer  :: HT_val(:)
    integer(C_INT),pointer  :: HT_outer(:),HT_inner(:)

contains
    !necessary t_H routines
    procedure :: set_from_Hcoo

    procedure :: add_child 
    procedure :: destroy_child    
    procedure :: copy_child 

    procedure :: optimize
    procedure :: mult_l


    procedure :: set_work_single
    procedure :: get_work_size_single

    procedure :: mult_r_disc

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
public t_H,t_H_eigen_mem
contains 

type(t_H_eigen_mem) function dummy_constructor()
    !might want some initialization for H and descr, but should work without
    !continue 
end function 

subroutine set_work_single(this,work,order)
    class(t_H_eigen_mem),intent(inout)      :: this
    class(work_ham_single),intent(inout)    :: work 
    integer,intent(in)                      :: order
    integer     :: sizes(2)
    integer     :: dim_mode

    Call this%t_H_eigen%set_work_single(work,order)
    if(this%row_max==0) ERROR STOP "cannot set work size of t_H_mkl_csr if row_max==0"
    Call this%get_work_size_single(sizes)
    Call work%set(sizes)
end subroutine

subroutine get_work_size_single(this,sizes)
    class(t_H_eigen_mem),intent(in) :: this
    integer,intent(out)             :: sizes(N_work_single)
    integer                         :: sizes_n(N_work_single)

    Call this%t_H_eigen%get_work_size_single(sizes_n)
    Call work_size_single(maxval(this%dim_l_single),this%row_max,sizes)
    sizes=max(sizes,sizes_n)
end subroutine

subroutine mult_r_disc(this,i_m,lat,N,ind_out,vec,ind_sum,ind_Mult,mat_mult,vec_mult)
    !Calculates the entries of the matrix * right vector product for the indices ind_out of the result vector
    USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
    ! input
    class(t_H_eigen_mem), intent(in)    :: this
    type(lattice), intent(in)           :: lat
    integer, intent(in)                 :: i_m          !index of the comp's right mode in the inner dim_mode
    integer, intent(in)                 :: N            !number of indices to calculated
    integer, intent(in)                 :: ind_out(N)   !indices to be calculated
    ! output
    real(8),intent(out)                 :: vec(N) ! dim_modes_inner(this%mode_r%order(comp))
    !temporary data
    integer,intent(inout)               :: ind_sum( N+1)                !ind_mult index where the a vec-entry start and end
    integer,intent(inout)               :: ind_mult(N*this%row_max)     !indices of the left array which have non-vanishing contributions to get the ind entries of the vec/mat product
    real(8),intent(inout)               :: mat_mult(N*this%row_max)     !matrix entries corresponding to ind_mult
    real(8),intent(inout)               :: vec_mult(N*this%row_max)     !values of discontiguous right mode which has to be evaluated (indices of ind_mult)

    !some local indices/ loop variables
    integer ::  i,j, i_outer,ii

    !get matrix indices and values whose vec.mat product constitute the output vector
    ind_sum(1)=0
    ii=0
    do i=1,N
        i_outer=ind_out(i)
        do j=this%HT_outer(i_outer)+1,this%HT_outer(i_outer+1)
            ii=ii+1
            mat_mult(ii)=this%HT_val(j)
            ind_mult(ii)=this%HT_inner(j)
        enddo
        ind_sum(i+1)=ii
    enddo
    ind_mult=ind_mult+1    !zero based index from c++

    !get the right vector values which are multiplied with the matrix
    Call this%mode_r%get_mode_disc(lat,ii,ind_mult(:ii),vec_mult(:ii))

    !multipy the matrix and right vector entries and sum together to the respective output entry
    vec_mult(:ii)=vec_mult(:ii)*mat_mult(:ii)
    do i=1,N
        vec(i)=sum(vec_mult(ind_sum(i)+1:ind_sum(i+1)))
    enddo
end subroutine 

!subroutine mult_r_ind(this,lat,N,ind_out,vec_out)
!    class(t_H_eigen_mem),intent(in) :: this
!    type(lattice),intent(in)        :: lat
!    integer,intent(in)              :: N
!    integer,intent(in)              :: ind_out(N)
!    real(8),intent(out)             :: vec_out(N)
!
!    real(8),pointer            :: modes(:)
!    real(8),allocatable,target :: vec(:)
!
!    Call this%mode_r%get_mode(lat,modes,vec)
!    Call eigen_H_mult_l_ind(this%HT,modes,N,ind_out,vec_out)
!end subroutine

subroutine mult_l(this,lat,res,work,alpha,beta)
    use m_derived_types, only: lattice
    class(t_H_eigen_mem),intent(in) :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)
    type(work_mode),intent(inout)   :: work
    real(8),intent(in),optional     :: alpha
    real(8),intent(in),optional     :: beta

    ! internal
    real(8),pointer            :: modes(:)
    real(8)                    :: alp, bet

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
    Call this%mode_l%get_mode(lat,modes,work)
    Call eigen_H_mult_r(this%HT,modes,res,alp,bet)
end subroutine 

subroutine optimize(this)
    class(t_H_eigen_mem),intent(inout)   :: this

    Call this%t_H_eigen%optimize()
    !THERE MIGHT BE SOMETHING TO DO HERE
    continue 
end subroutine

subroutine copy_child(this,Hout)
    class(t_H_eigen_mem),intent(in)     :: this
    class(t_H_base),intent(inout)        :: Hout
    
    Call this%t_H_eigen%copy_child(Hout)
    select type(Hout)
    class is(t_H_eigen_mem)
        Call eigen_H_copy(this%HT,Hout%HT) 
        Call Hout%set_auxiliaries()
    class default
        STOP "Cannot copy t_H_eigen_mem type with Hamiltonian that is not a class of t_H_eigen_mem"
    end select
end subroutine

subroutine add_child(this,H_in)
    class(t_H_eigen_mem),intent(inout)    :: this
    class(t_H_base),intent(in)            :: H_in

    Call this%t_H_eigen%add_child(H_in)

    select type(H_in)
    class is(t_H_eigen_mem)
        Call eigen_H_add(this%HT,H_in%HT)
        Call this%set_auxiliaries()
    class default
        STOP "Cannot add t_H_eigen_mem type with Hamiltonian that is not a class of t_H_eigen_mem"
    end select
end subroutine 

subroutine destroy_child(this)
    class(t_H_eigen_mem),intent(inout)    :: this

    Call this%t_H_eigen%destroy_child()
    if(this%is_set())then
        Call eigen_H_destroy(this%HT)
        nullify(this%HT_val,this%HT_outer,this%HT_inner)
        this%row_max=0
        this%dim_l_single=0
    endif
end subroutine

subroutine set_from_Hcoo(this,H_coo)
    class(t_H_eigen_mem),intent(inout)  :: this
    type(t_H_coo),intent(inout)     :: H_coo

    Call this%t_H_eigen%set_from_Hcoo(H_coo)
    Call eigen_get_transpose(this%H,this%HT)
    Call this%set_auxiliaries()
end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!            MPI ROUTINES           !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine send(this,ithread,tag,com)
    use mpi_basic                
    class(t_H_eigen_mem),intent(in) :: this
    integer,intent(in)              :: ithread
    integer,intent(in)              :: tag
    integer,intent(in)              :: com

#ifdef CPP_MPI
    Call this%t_H_eigen%send(ithread,tag,com)
    Call eigen_H_send(ithread,tag,this%HT,com) 
#else
    continue
#endif
end subroutine

subroutine recv(this,ithread,tag,com)
    use mpi_basic                
    class(t_H_eigen_mem),intent(inout)  :: this
    integer,intent(in)                  :: ithread
    integer,intent(in)                  :: tag
    integer,intent(in)                  :: com

#ifdef CPP_MPI
    Call this%t_H_eigen%recv(ithread,tag,com)
    Call eigen_H_recv(ithread,tag,this%HT,com) 
    Call this%set_auxiliaries()
#else
    continue
#endif
end subroutine

subroutine bcast(this,comm)
    use mpi_basic                
    class(t_H_eigen_mem),intent(inout)  ::  this
    type(mpi_type),intent(in)           ::  comm
#ifdef CPP_MPI
    Call this%t_H_eigen%bcast(comm)
    Call eigen_H_bcast(comm%id,comm%mas,comm%ismas,this%HT,comm%com) 
    if(.not.comm%ismas) Call this%set_auxiliaries()
#else
    continue
#endif
end subroutine 

subroutine distribute(this,comm)
    use mpi_basic                
    class(t_H_eigen_mem),intent(inout)  ::  this
    type(mpi_type),intent(in)           ::  comm
#ifdef CPP_MPI
    Call this%t_H_eigen%distribute(comm)
    Call eigen_H_distribute(comm%id,comm%mas,comm%ismas,this%HT,comm%com)
    if(.not.comm%ismas) Call this%set_auxiliaries()
#else
    continue
#endif
end subroutine 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           HELPER FUNCTIONS                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_auxiliaries(this)
    class(t_H_eigen_mem),intent(inout)    :: this

    Call this%t_H_eigen%set_auxiliaries()
    Call set_H_ptr(this)
    Call set_row_max(this)
end subroutine

subroutine set_H_ptr(this)
    class(t_H_eigen_mem),intent(inout)    :: this

    type(C_PTR)     :: col,row,val
    integer         :: dimH(2), nnz

    nullify(this%HT_inner,this%HT_outer,this%HT_val)
    Call eigen_get_dat(this%HT,nnz,dimH,col,row,val)
    if(nnz/=this%nnz) ERROR STOP "number of entries of Hamiltonian and Hamiltonian transpose should be identical"
    CAll C_F_POINTER(col, this%HT_outer, [dimH(2)+1])
    CAll C_F_POINTER(row, this%HT_inner, [nnz]) 
    CAll C_F_POINTER(val, this%HT_val, [nnz]) 
end subroutine

subroutine set_row_max(this)
    class(t_H_eigen_mem),intent(inout)  :: this
    !variable to get multiplication for single evaluation (estimate sizes...)
    integer,allocatable                 :: tmp(:)
    integer                             :: i

    allocate(tmp(size(this%HT_outer)-1))
    do i=1,size(tmp)
        tmp(i)=this%HT_outer(i+1)-this%HT_outer(i)
    enddo
    this%row_max=maxval(tmp)
end subroutine
#endif 
end module
