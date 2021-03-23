module m_H_eigen
#ifdef CPP_EIGEN_H
!Hamiltonian type specifications using dense matrices and no external library
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
    procedure :: mult_r,mult_l
    !overriding mult_l/r_single for better efficiency
    procedure :: mult_r_single,mult_l_single
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

subroutine mult_r(this,lat,res)
    !mult
    use m_derived_types, only: lattice
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
    class(t_H),intent(inout)        :: Hout
    
    select type(Hout)
    class is(t_H_eigen)
        Call eigen_H_copy(this%H,Hout%H) 
    class default
        STOP "Cannot copy t_H_eigen type with Hamiltonian that is not a class of t_H_eigen"
    end select
end subroutine

subroutine add_child(this,H_in)
    class(t_H_eigen),intent(inout)    :: this
    class(t_H),intent(in)             :: H_in

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
    ! internal
    real(8),pointer                 :: modes_l(:),modes_r(:)
    real(8),allocatable,target      :: vec_l(:),vec_r(:)
    integer                         :: bnd(2)
    integer                         :: size_vec_r

    Call this%mode_l%get_mode(lat,modes_l,vec_l)
    ERROR STOP "THIS PROBABLY NO LONGER WORKS WITH THE NEW MODE_L/MODE_R"   !and in general might be much more difficult to implement with eg. rank 4 in M-space only
!    Call lat%point_order(this%op_l,this%dimH(1),modes_l,vec_l)
    Call lat%point_order_single(this%op_r,i_m,dim_bnd,this%dim_mode(2),modes_r,vec_r,bnd)

    size_vec_r=bnd(2)-bnd(1)+1
    Call eigen_H_eval_single(bnd(1)-1,size_vec_r,modes_l,modes_r,this%H,E)
    
    if(allocated(vec_l)) deallocate(vec_l)
    if(allocated(vec_r)) deallocate(vec_r)
end subroutine 

#endif 
end module
