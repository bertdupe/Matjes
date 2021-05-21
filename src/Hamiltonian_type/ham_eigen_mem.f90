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
    integer(C_INT),pointer  :: HT_col(:),HT_row(:)
    !helper variables
    integer                 :: row_max=0 !maximal number of entries per col

contains
    !necessary t_H routines
    procedure :: set_from_Hcoo

    procedure :: add_child 
    procedure :: destroy_child    
    procedure :: copy_child 

    procedure :: optimize
    procedure :: mult_l
!    procedure :: mult_l_cont
!    procedure :: mult_l_disc
!    procedure :: mult_r_single,mult_l_single   !can probably be made more efficient, but I am not sure if it is still necessary
    procedure :: mult_r_ind 
!    procedure :: mult_r_disc_disc
!    procedure :: get_ind_mult_l


    procedure :: set_work_single
    procedure :: get_work_size_single

    procedure :: mult_r_single

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
    if(.not.any(order==this%op_l))then
        if(any(order==this%op_r))then
            ERROR STOP "DECIDE HOW TO TREATE CASE WITHOUT ENTRY"
        else
            ERROR STOP "Hamiltonian has no component of considered single energy evaluation, take it out or consider it somehow else"
        endif
    endif
    Call this%get_work_size_single(sizes)
    Call work%set(sizes)
end subroutine

subroutine get_work_size_single(this,sizes)
    class(t_H_eigen_mem),intent(in) :: this
    integer,intent(out)             :: sizes(N_work_single)
    integer                         :: sizes_n(N_work_single)

    Call this%t_H_eigen%get_work_size_single(sizes_n)
    Call work_size_single(this%dim_l_single,this%row_max,sizes)
    sizes=max(sizes,sizes_n)
end subroutine

subroutine mult_r_single(this,i_m,comp,lat,work,vec)
    !Calculates the entries of the matrix * right vector product which corresponds to the i_m's site of component comp of the left modes
    USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
    ! input
    class(t_H_eigen_mem), intent(in)    :: this
    type(lattice), intent(in)           :: lat
    integer, intent(in)                 :: i_m           !index of the comp's right mode in the inner dim_mode
    integer, intent(in)                 :: comp          !component of right mode
    !temporary data
    type(work_ham_single),intent(inout) :: work          !data type containing the temporary data for this calculation to prevent constant allocations/deallocations
    ! output
    real(8),intent(inout)               :: vec(:)        !dim_modes_inner(this%mode_l%order(comp))

    integer :: i,j,i_col
    integer :: ii

#ifdef CPP_USE_WORK
    !temporary data slices
    integer,pointer,contiguous          :: ind_out(:)    !indices of all right mode entries which contain the order paramtere order of the site corresponding to i_m
    integer,pointer,contiguous          :: ind_sum(:)    !ind_mult index where the a vec-entry start and end
    integer,pointer,contiguous          :: ind_mult(:)   !indices of the left array which have non-vanishing contributions to get the ind entries of the vec/mat product
    real(8),pointer,contiguous          :: mat_mult(:)   !matrix entries corresponding to ind_mult
    real(8),pointer,contiguous          :: vec_r(:)      !values of discontiguous right mode which has to be evaluated (indices of ind_mult)


    !associate temporary arrays
    !!int vector slices
    ind_out (1:this%dim_r_single             )=>work%int_arr (1                               :this%dim_l_single                    )
    ind_sum (1:this%dim_r_single+1           )=>work%int_arr (1+this%dim_r_single             :this%dim_r_single* 2               +1)
    ind_mult(1:this%dim_r_single*this%col_max)=>work%int_arr (1+this%dim_r_single*2+1         :this%dim_r_single*(2+ this%col_max)+1)
    !!real vector slices
    vec_r   (1:this%dim_r_single*this%col_max)=>work%real_arr(1                               :this%dim_r_single*this%col_max                  )
    mat_mult(1:this%dim_r_single*this%col_max)=>work%real_arr(1+this%dim_r_single*this%col_max:this%dim_r_single*this%col_max*2                )
#else
    !temporary arrays
    integer                     :: ind_out (this%dim_l_single             )     !indices of all right mode entries which contain the order paramtere order of the site corresponding to i_m
    integer                     :: ind_sum (this%dim_l_single+1           )     !ind_mult index where the a vec-entry start and end
    integer                     :: ind_mult(this%dim_l_single*this%row_max)     !indices of the left array which have non-vanishing contributions to get the ind entries of the vec/mat product
    real(8)                     :: mat_mult(this%dim_l_single*this%row_max)     !matrix entries corresponding to ind_mult
    real(8)                     :: vec_r   (this%dim_l_single*this%row_max)     !values of discontiguous right mode which has to be evaluated (indices of ind_mult)
#endif
    !get indices of the output vector 
    Call this%mode_l%get_ind_site_expl(comp,i_m,this%dim_l_single,ind_out)

    !get matrix indices and values whose vec.mat product constitute the output vector
    ind_sum(1)=0
    ii=0
    do i=1,this%dim_l_single
        i_col=ind_out(i)
        do j=this%HT_col(i_col)+1,this%HT_col(i_col+1)
            ii=ii+1
            mat_mult(ii)=this%HT_val(j)
            ind_mult(ii)=this%HT_row(j)
        enddo
        ind_sum(i+1)=ii
    enddo
    ind_mult=ind_mult+1    !zero based index from c++

    !get the right vector values which are multiplied with the matrix
    Call this%mode_r%get_mode_disc_expl(lat,ii,ind_mult(:ii),vec_r(:ii))

    !multipy the matrix and right vector entries and sum together to the respective output entry
    vec_r(:ii)=vec_r(:ii)*mat_mult(:ii)
    do i=1,this%dim_l_single
        vec(i)=sum(vec_r(ind_sum(i)+1:ind_sum(i+1)))
    enddo
end subroutine 

!subroutine get_ind_mult_l(this,ind_in,N_out,ind_out)
!    !get the indicies in ind_out(:N_out) of the right vector which are necessary to 
!    !get all components that 
!    class(t_H_eigen_mem),intent(in) :: this
!    integer,intent(in)              :: ind_in(:)
!    integer,intent(out)             :: N_out
!    integer,intent(inout)           :: ind_out(:)
!
!    N_out=size(ind_out)
!    Call eigen_H_get_ind_mult_r(this%HT,size(ind_in),ind_in,N_out,ind_out)
!    ind_out(1:N_out)=ind_out(1:N_out)+1
!end subroutine 

!subroutine mult_l_cont(this,bnd,vec,res)
!    !multiply to the right with a continuous section of the right vector
!    class(t_H_eigen_mem),intent(in) :: this
!    integer,intent(in)              :: bnd(2)
!    real(8),intent(in)              :: vec(bnd(2)-bnd(1)+1)
!    real(8),intent(inout)           :: res(:)   !result matrix-vector product
!
!    Call eigen_H_mult_mat_vec_cont(this%HT,bnd(1),bnd(2),vec,res)
!end subroutine 
!
!subroutine mult_l_disc(this,N,ind,vec,res)
!    !multiply to the right with a discontinuous section of the right vector
!    class(t_H_eigen_mem),intent(in) :: this
!    integer,intent(in)              :: N
!    integer,intent(in)              :: ind(N)
!    real(8),intent(in)              :: vec(N)
!    real(8),intent(inout)           :: res(:)   !result matrix-vector product
!
!    Call eigen_H_mult_mat_vec_disc(this%HT,N,ind,vec,res)
!end subroutine 

subroutine mult_r_ind(this,lat,N,ind_out,vec_out)
    class(t_H_eigen_mem),intent(in) :: this
    type(lattice),intent(in)        :: lat
    integer,intent(in)              :: N
    integer,intent(in)              :: ind_out(N)
    real(8),intent(out)             :: vec_out(N)

    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)

    Call this%mode_r%get_mode(lat,modes,vec)
    Call eigen_H_mult_l_ind(this%HT,modes,N,ind_out,vec_out)
end subroutine

subroutine mult_l(this,lat,res)
    use m_derived_types, only: lattice
    class(t_H_eigen_mem),intent(in) :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)
    ! internal
    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)

    Call this%mode_l%get_mode(lat,modes,vec)
    Call eigen_H_mult_mat_vec(this%HT,modes,res)
    if(allocated(vec)) deallocate(vec)
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
        nullify(this%HT_val,this%HT_col,this%HT_row)
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

    nullify(this%HT_row,this%HT_col,this%HT_val)
    Call eigen_get_dat(this%HT,nnz,dimH,col,row,val)
    if(nnz/=this%nnz) ERROR STOP "number of entries of Hamiltonian and Hamiltonian transpose should be identical"
    CAll C_F_POINTER(col, this%HT_col, [dimH(2)+1])
    CAll C_F_POINTER(row, this%HT_row, [nnz]) 
    CAll C_F_POINTER(val, this%HT_val, [nnz]) 
end subroutine

subroutine set_row_max(this)
    class(t_H_eigen_mem),intent(inout)  :: this
    !variable to get multiplication for single evaluation (estimate sizes...)
    integer,allocatable                 :: tmp(:)
    integer                             :: i

    allocate(tmp(size(this%HT_col)-1))
    do i=1,size(tmp)
        tmp(i)=this%HT_col(i+1)-this%HT_col(i)
    enddo
    this%row_max=maxval(tmp)
end subroutine
#endif 
end module
