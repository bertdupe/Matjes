module m_H_eigen
#ifdef CPP_EIGEN_H
!Hamiltonian type specifications using Eigen to save and evaluate the Hamiltonian
use m_derived_types, only: lattice, number_different_order_parameters
use m_eigen_H_interface
use m_H_coo_based
use m_mode_construction_rankN_sparse_col

type,extends(t_H_coo_based) :: t_H_eigen
    type(C_PTR)     ::  H=c_null_ptr

!some pointer to Hamiltonian
contains
    !necessary t_H routines
    procedure :: eval_single

    procedure :: set_from_Hcoo

    procedure :: add_child 
    procedure :: bcast_child 
    procedure :: destroy_child    
    procedure :: copy_child 
    procedure :: optimize

    procedure :: mult_l,           mult_r
    procedure :: mult_l_cont,      mult_r_cont
    procedure :: mult_l_disc,      mult_r_disc
    procedure :: mult_r_single,    mult_l_single
    procedure :: mult_l_ind,       mult_r_ind
    procedure :: mult_l_disc_disc, mult_r_disc_disc
    procedure :: get_ind_mult_l,   get_ind_mult_r
end type

interface t_H_mkl_csr
    procedure :: dummy_constructor
end interface 

private
public t_H,t_H_eigen
contains 

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
    class(t_H_base),intent(inout)        :: Hout
    
    select type(Hout)
    class is(t_H_eigen)
        Call eigen_H_copy(this%H,Hout%H) 
        Call this%copy_deriv(Hout)
    class default
        STOP "Cannot copy t_H_eigen type with Hamiltonian that is not a class of t_H_eigen"
    end select
end subroutine

subroutine add_child(this,H_in)
    class(t_H_eigen),intent(inout)    :: this
    class(t_H_base),intent(in)             :: H_in

    select type(H_in)
    class is(t_H_eigen)
        Call eigen_H_add(this%H,H_in%H)
    class default
        STOP "Cannot add t_H_eigen type with Hamiltonian that is not a class of t_H_eigen"
    end select
end subroutine 

subroutine bcast_child(this,comm)
    use mpi_basic                
    class(t_H_eigen),intent(inout)  ::  this
    type(mpi_type),intent(in)       ::  comm
#ifdef CPP_MPI
    Call eigen_H_bcast(comm%id,comm%mas,comm%ismas,this%H,comm%com) 
#else
    continue
#endif
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
end subroutine 

subroutine eval_single(this,E,i_m,dim_bnd,lat)
    use m_derived_types, only: lattice
    ! input
    class(t_H_eigen),intent(in)     :: this
    type(lattice), intent(in)       :: lat
    integer, intent(in)             :: i_m
    integer, intent(in)             :: dim_bnd(2,number_different_order_parameters)
    ! output
    real(8), intent(out)            :: E
    integer,allocatable             :: ind(:)
    real(8),allocatable             :: vec(:)

    real(8),allocatable             :: vec_out(:)
    integer,allocatable             :: ind_out(:)
    real(8),allocatable             :: vec_l(:)
    integer                         :: N_out

    !dim_order_bnd...
    Call this%mode_r%get_mode_single_disc(lat,1,i_m,ind,vec)
    N_out=size(ind)*10  !arbitrary, hopefully large enough (otherwise eigen_H_milt_mat_disc_disc crashes
    allocate(vec_out(N_out), ind_out(N_out))
    Call eigen_H_mult_mat_disc_disc(this%H,size(ind),ind,vec,N_out,ind_out,vec_out)
    Call this%mode_l%get_mode_disc(lat,ind_out(:N_out),vec_l)
    E=DOT_PRODUCT(vec_l(:N_out),vec_out(:N_out))
end subroutine 
#endif 
end module
