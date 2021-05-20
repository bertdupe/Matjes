module m_H_eigen
#ifdef CPP_EIGEN
!Hamiltonian type specifications using Eigen to save and evaluate the Hamiltonian
use m_derived_types, only: lattice, number_different_order_parameters
use m_eigen_H_interface
use m_H_coo_based
use m_mode_construction_rankN_sparse_col
use m_work_ham_single

type,extends(t_H_coo_based) :: t_H_eigen
    type(C_PTR)     ::  H=c_null_ptr  !pointer to sparse matrix c-object 

    !pointers to Hamiltonian data handled by eigen (col major format) (zero based)
    real(C_DOUBLE),pointer  :: mat_val(:)
    integer(C_INT),pointer  :: mat_col(:),mat_row(:)
    integer                 :: nnz
    !helper variables
    integer                 :: col_max=0 !maximal number of entries per col
    integer                 :: dim_single !dimension for inner work size array of set order (set_work_size_single)
contains
    !necessary t_H routines
    procedure :: eval_single

    procedure :: set_from_Hcoo

    procedure :: add_child 
    procedure :: destroy_child
    procedure :: copy_child 
    procedure :: optimize

    procedure :: mult_l,           mult_r
    procedure :: mult_l_cont,      mult_r_cont
    procedure :: mult_l_disc,      mult_r_disc
    procedure :: mult_l_single,    mult_r_single
    procedure :: mult_l_ind,       mult_r_ind
    procedure :: mult_l_disc_disc, mult_r_disc_disc
    procedure :: get_ind_mult_l,   get_ind_mult_r

    procedure :: set_work_size_single

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

subroutine set_work_size_single(this,work,order)
    class(t_H_eigen),intent(inout)          :: this
    class(work_ham_single),intent(inout)    :: work 
    integer,intent(in)                      :: order
    integer     :: sizes(2)
    integer     :: dim_mode

    if(.not.this%is_set()) ERROR STOP "cannot set work size of hamiltonian if it is not set"
    if(this%col_max==0) ERROR STOP "cannot set work size of t_H_mkl_csr if row_max==0"
    if(.not.any(order==this%op_r))then
        if(any(order==this%op_l))then
            ERROR STOP "So far cannot consider Hamiltonian which has the single evaluation operator in the left side, but not on the right (some csr implementation?)"
        else
            ERROR STOP "Hamiltonian has no component of considered single energy evaluation, take it out or consider it somehow else"
        endif
    endif
    Call this%mode_r%get_mode_single_size(order,this%dim_single)    !multiply to the right for single evaluation as that is easier in csc format
    sizes(1)= this%dim_single*(2*this%col_max+1)
    sizes(2)= this%dim_single*(  this%col_max+1)
    Call work%set(sizes)
end subroutine

subroutine eval_single(this,E,i_m,order,lat,work)
    USE, INTRINSIC :: ISO_C_BINDING , ONLY : C_INT, C_DOUBLE
    ! input
    class(t_H_eigen), intent(in)        :: this
    type(lattice), intent(in)           :: lat
    integer, intent(in)                 :: i_m
    integer, intent(in)                 :: order
    ! output
    real(8), intent(out)                :: E
    !temporary data
    type(work_ham_single),intent(inout) :: work    !data type containing the temporary data for this calculation to prevent constant allocations/deallocations
    !temporary data slices
    integer,pointer,contiguous          :: ind(:)       !indices of all right mode entries which contain the order parameter order of the site corresponding to i_m
    real(8),pointer,contiguous          :: vec(:)       !values corresponding to ind
    integer,pointer,contiguous          :: ind_out(:)   !indices of the result array multipling the vector (vec/ind) to the matrix
    real(8),pointer,contiguous          :: vec_out(:)   !values corresponding to ind_out
    real(8),pointer,contiguous          :: vec_l(:)     !values of discontiguous mode array on right side (indices of ind_out)

    !some local indices/ loop variables
    integer ::  i,j, i_col
    integer :: ii

    !associate temporary arrays
    ind    (1:this%dim_single             )=>work%int_arr (1                                 :this%dim_single                   )
    ind_out(1:this%dim_single*this%col_max)=>work%int_arr (1+this%dim_single                 :this%dim_single*(1+  this%col_max))
    vec    (1:this%dim_single             )=>work%real_arr(1                                 :this%dim_single                   )
    vec_out(1:this%dim_single*this%col_max)=>work%real_arr(1+this%dim_single                 :this%dim_single*(1+  this%col_max))
    vec_l  (1:this%dim_single*this%col_max)=>work%real_arr(1+this%dim_single*(1+this%col_max):this%dim_single*(1+2*this%col_max))

    !get right mode corresponding to site i_m of order order
    Call this%mode_r%get_mode_single(lat,1,i_m,this%dim_single,ind,vec)    !get this to work with different orders (1 is not order here but component of left mode)

    !Calculate matrix,right mode product only for the necessary discontiguous mode-indices
    ii=0
    do i=1,size(ind)
        i_col=ind(i)
        do j=this%mat_col(i_col)+1,this%mat_col(i_col+1)
            ii=ii+1
            vec_out(ii)=this%mat_val(j)*vec(i)
            ind_out(ii)=this%mat_row(j)
        enddo
    enddo
    ind_out=ind_out+1   !zero based index from c++
    
    !get left mode for indices of the mat/vec product 
    Call this%mode_r%get_mode_disc_expl(lat,ii,ind_out(:ii),vec_l(:ii))

    !Get the energy
    E=DOT_PRODUCT(vec_l(:ii),vec_out(:ii))
    nullify(ind,ind_out,vec,vec_out,vec_l)
end subroutine 


type(t_H_eigen) function dummy_constructor()
    !might want some initialization for H and descr, but should work without
    !continue 
end function 

subroutine mult_l_disc_disc(this,ind_l,vec_l,ind_r,vec_out)
    class(t_H_eigen),intent(in)     :: this
    integer,intent(in)              :: ind_l(:)
    real(8),intent(in)              :: vec_l(:)
    integer,intent(in)              :: ind_r(:)
    real(8),intent(inout)           :: vec_out(:)

    Call eigen_mult_l_disc_disc(this%H,size(ind_l),ind_l,vec_l,size(ind_r),ind_r,vec_out)
end subroutine

subroutine mult_r_disc_disc(this,ind_r,vec_r,ind_l,vec_out)
    class(t_H_eigen),intent(in)     :: this
    integer,intent(in)              :: ind_r(:)
    real(8),intent(in)              :: vec_r(:)
    integer,intent(in)              :: ind_l(:)
    real(8),intent(inout)           :: vec_out(:)

    Call eigen_mult_r_disc_disc(this%H,size(ind_r),ind_r,vec_r,size(ind_l),ind_l,vec_out)
end subroutine

subroutine get_ind_mult_r(this,ind_in,N_out,ind_out)
    !get the indicies in ind_out(:N_out) of the right vector which are necessary to 
    !get all components that 
    class(t_H_eigen),intent(in)     :: this
    integer,intent(in)              :: ind_in(:)
    integer,intent(out)             :: N_out
    integer,intent(inout)           :: ind_out(:)

    N_out=size(ind_out)
    Call eigen_H_get_ind_mult_r(this%H,size(ind_in),ind_in,N_out,ind_out)
    ind_out(1:N_out)=ind_out(1:N_out)+1
end subroutine 

subroutine get_ind_mult_l(this,ind_in,N_out,ind_out)
    !get the indicies in ind_out(:N_out) of the right vector which are necessary to 
    !get all components that 
    class(t_H_eigen),intent(in)     :: this
    integer,intent(in)              :: ind_in(:)
    integer,intent(out)             :: N_out
    integer,intent(inout)           :: ind_out(:)

    N_out=size(ind_out)
    Call eigen_H_get_ind_mult_l(this%H,size(ind_in),ind_in,N_out,ind_out)
    ind_out(1:N_out)=ind_out(1:N_out)+1
end subroutine 

subroutine mult_r(this,lat,res)
    !mult
    class(t_H_eigen),intent(in)     :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)   !result matrix-vector product
    ! internal
    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)

    Call this%mode_r%get_mode(lat,modes,vec)
    Call eigen_H_mult_mat_vec(this%H,modes,res)
    if(allocated(vec)) deallocate(vec)
end subroutine 

subroutine mult_l_cont(this,bnd,vec,res)
    !multiply to the right with a continuous section of the right vector
    class(t_H_eigen),intent(in)     :: this
    integer,intent(in)              :: bnd(2)
    real(8),intent(in)              :: vec(bnd(2)-bnd(1)+1)
    real(8),intent(inout)           :: res(:)   !result matrix-vector product

    Call eigen_H_mult_vec_mat_cont(this%H,bnd(1),bnd(2),vec,res)
end subroutine 

subroutine mult_r_cont(this,bnd,vec,res)
    !multiply to the right with a continuous section of the right vector
    class(t_H_eigen),intent(in)     :: this
    integer,intent(in)              :: bnd(2)
    real(8),intent(in)              :: vec(bnd(2)-bnd(1)+1)
    real(8),intent(inout)           :: res(:)   !result matrix-vector product

    Call eigen_H_mult_mat_vec_cont(this%H,bnd(1),bnd(2),vec,res)
end subroutine 

subroutine mult_l_disc(this,N,ind,vec,res)
    !multiply to the right with a discontinuous section of the right vector
    class(t_H_eigen),intent(in)     :: this
    integer,intent(in)              :: N
    integer,intent(in)              :: ind(N)
    real(8),intent(in)              :: vec(N)
    real(8),intent(inout)           :: res(:)   !result matrix-vector product

    Call eigen_H_mult_vec_mat_disc(this%H,N,ind,vec,res)
end subroutine 

subroutine mult_r_disc(this,N,ind,vec,res)
    !multiply to the right with a discontinuous section of the right vector
    class(t_H_eigen),intent(in)     :: this
    integer,intent(in)              :: N
    integer,intent(in)              :: ind(N)
    real(8),intent(in)              :: vec(N)
    real(8),intent(inout)           :: res(:)   !result matrix-vector product

    Call eigen_H_mult_mat_vec_disc(this%H,N,ind,vec,res)
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


subroutine mult_r_single(this,i_site,lat,res)
    !mult
    use m_derived_types, only: lattice
    class(t_H_eigen),intent(in)     :: this
    integer,intent(in)              :: i_site
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)   !result matrix-vector product
    ! internal
    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)
    integer                    :: bnd(2)

    Call this%mode_r%get_mode(lat,modes,vec)
    bnd(1)=this%dim_mode(1)*(i_site-1)
    bnd(2)=this%dim_mode(1)*(i_site)-1
    Call eigen_H_mult_mat_vec_single(this%H,bnd(1),bnd(2),modes,res)
    if(allocated(vec)) deallocate(vec)
end subroutine 


subroutine mult_l_single(this,i_site,lat,res)
    !mult
    use m_derived_types, only: lattice
    class(t_H_eigen),intent(in)     :: this
    integer,intent(in)              :: i_site
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)   !result matrix-vector product
    ! internal
    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)
    integer                    :: bnd(2)

    Call this%mode_l%get_mode(lat,modes,vec)
    bnd(1)=this%dim_mode(2)*(i_site-1)
    bnd(2)=this%dim_mode(2)*(i_site)-1
    Call eigen_H_mult_vec_mat_single(this%H,bnd(1),bnd(2),modes,res)
    if(allocated(vec)) deallocate(vec)
end subroutine 


subroutine mult_l(this,lat,res)
    use m_derived_types, only: lattice
    class(t_H_eigen),intent(in)     :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)
    ! internal
    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)

    Call this%mode_l%get_mode(lat,modes,vec)
    Call eigen_H_mult_vec_mat(this%H,modes,res)
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
        Call set_auxiliaries(this)
    class default
        STOP "Cannot add t_H_eigen type with Hamiltonian that is not a class of t_H_eigen"
    end select
end subroutine 

subroutine destroy_child(this)
    class(t_H_eigen),intent(inout)    :: this

    if(this%is_set())then
        Call eigen_H_destroy(this%H)
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
    Call set_auxiliaries(this)
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
    Call set_auxiliaries(this)
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
    if(.not.comm%ismas) Call set_auxiliaries(this)
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
    if(.not.comm%ismas) Call set_auxiliaries(this)
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

    nullify(this%mat_row,this%mat_col,this%mat_val)
    Call eigen_get_dat(this%H,this%nnz,dimH,col,row,val)
    CAll C_F_POINTER(col, this%mat_col, [dimH(1)+1])
    CAll C_F_POINTER(row, this%mat_row, [this%nnz]) 
    CAll C_F_POINTER(val, this%mat_val, [this%nnz]) 
end subroutine

subroutine set_col_max(this)
    class(t_H_eigen),intent(inout)    :: this
    !variable to get multiplication for single evaluation (estimate sizes...)
    integer,allocatable               :: tmp(:)
    integer                           :: i

    allocate(tmp(size(this%mat_col)-1))
    do i=1,size(tmp)
        tmp(i)=this%mat_col(i+1)-this%mat_col(i)
    enddo
    this%col_max=maxval(tmp)
end subroutine
#endif 
end module
