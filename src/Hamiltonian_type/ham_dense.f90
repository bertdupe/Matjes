module m_H_dense
!Hamiltonian type specifications using dense matrices and no external library
use m_derived_types, only: lattice, number_different_order_parameters
use m_H_coo_based

type,extends(t_H_coo_based) :: t_H_dense
    real(8),allocatable   :: H(:,:)
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
    procedure :: mult_l_cont,mult_r_cont
    procedure :: mult_l_disc,mult_r_disc
    procedure :: mult_r_single,mult_l_single
end type

private
public t_H,t_H_dense
contains 

subroutine mult_r(this,lat,res)
    !mult
    use m_derived_types, only: lattice
    class(t_h_dense),intent(in)   :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)   !result matrix-vector product
    ! internal
    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)

    Call this%mode_r%get_mode(lat,modes,vec)
    if(size(res)/=this%dimH(1)) STOP "size of vec is wrong"
    res=matmul(this%H,modes)
    if(allocated(vec)) deallocate(vec)
end subroutine 

subroutine mult_l(this,lat,res)
    use m_derived_types, only: lattice
    class(t_h_dense),intent(in)   :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)
    ! internal
    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)

    Call this%mode_l%get_mode(lat,modes,vec)
    if(size(res)/=this%dimH(2)) STOP "size of vec is wrong"
    res=matmul(modes,this%H)
    if(allocated(vec)) deallocate(vec)
end subroutine 

subroutine mult_r_cont(this,bnd,vec,res)
    class(t_H_dense),intent(in)    :: this
    integer,intent(in)           :: bnd(2)
    real(8),intent(in)           :: vec(bnd(2)-bnd(1)+1)
    real(8),intent(inout)        :: res(:)   !result matrix-vector product

    STOP "IMPLEMENT if necessary"
end subroutine 

subroutine mult_l_cont(this,bnd,vec,res)
    class(t_H_dense),intent(in)    :: this
    integer,intent(in)           :: bnd(2)
    real(8),intent(in)           :: vec(bnd(2)-bnd(1)+1)
    real(8),intent(inout)        :: res(:)   !result matrix-vector product

    STOP "IMPLEMENT if necessary"
end subroutine 

subroutine mult_r_disc(this,N,ind,vec,res)
    class(t_H_dense),intent(in)    :: this
    integer,intent(in)           :: N
    integer,intent(in)           :: ind(N)
    real(8),intent(in)           :: vec(N)
    real(8),intent(inout)        :: res(:)   !result matrix-vector product

    STOP "IMPLEMENT if necessary"
end subroutine 

subroutine mult_l_disc(this,N,ind,vec,res)
    class(t_H_dense),intent(in)    :: this
    integer,intent(in)           :: N
    integer,intent(in)           :: ind(N)
    real(8),intent(in)           :: vec(N)
    real(8),intent(inout)        :: res(:)   !result matrix-vector product

    STOP "IMPLEMENT if necessary"
end subroutine 



subroutine mult_r_single(this,i_site,lat,res)
    !mult
    use m_derived_types, only: lattice
    class(t_h_dense),intent(in) :: this
    integer,intent(in)          :: i_site
    type(lattice), intent(in)   :: lat
    real(8), intent(inout)      :: res(:)   !result matrix-vector product
    ! internal
    integer                     :: bnd(2)
    real(8),pointer             :: modes(:)
    real(8),allocatable,target  :: vec(:)

    Call this%mode_r%get_mode(lat,modes,vec)
    if(size(res)/=this%dimH(1)) STOP "size of vec is wrong"
    bnd(1)=this%dim_mode(1)*(i_site-1)+1
    bnd(2)=this%dim_mode(1)*(i_site)
    !terrible implementation striding-wise
    res=matmul(this%H(bnd(1):bnd(2),:),modes)

    if(allocated(vec)) deallocate(vec)
end subroutine 

subroutine mult_l_single(this,i_site,lat,res)
    !mult
    use m_derived_types, only: lattice
    class(t_h_dense),intent(in) :: this
    integer,intent(in)          :: i_site
    type(lattice), intent(in)   :: lat
    real(8), intent(inout)      :: res(:)   !result matrix-vector product
    ! internal
    integer                     :: bnd(2)
    real(8),pointer             :: modes(:)
    real(8),allocatable,target  :: vec(:)

    Call this%mode_l%get_mode(lat,modes,vec)
    if(size(res)/=this%dimH(1)) STOP "size of vec is wrong"
    bnd(1)=this%dim_mode(2)*(i_site-1)+1
    bnd(2)=this%dim_mode(2)*(i_site)
    res=matmul(modes,this%H(:,bnd(1):bnd(2)))

    if(allocated(vec)) deallocate(vec)
end subroutine 

subroutine optimize(this)
    class(t_h_dense),intent(inout)   :: this

    !nothing to optimize here
    continue 
end subroutine


subroutine copy_child(this,Hout)
    class(t_h_dense),intent(in)   :: this
    class(t_H_base),intent(inout)        :: Hout
    
    select type(Hout)
    class is(t_h_dense)
        allocate(Hout%H,source=this%H)
        Call this%copy_deriv(Hout)
    class default
        STOP "Cannot copy t_h_dense type with Hamiltonian that is not a class of t_h_dense"
    end select
end subroutine

subroutine bcast_child(this,comm)
    use mpi_basic                
    class(t_H_dense),intent(inout)        ::  this
    type(mpi_type),intent(in)       ::  comm
#ifdef CPP_MPI
    integer     :: ierr
    integer     :: N(2)
    
    if(comm%ismas)then
        N=shape(this%H)
    endif
    Call MPI_Bcast(N, 2, MPI_INTEGER, comm%mas, comm%com,ierr)
    if(.not.comm%ismas)then
        allocate(this%H(N(1),N(2)))
    endif
    Call MPI_Bcast(this%H,size(this%H), MPI_REAL8, comm%mas, comm%com,ierr)
#else
    continue
#endif
end subroutine 

subroutine add_child(this,H_in)
    class(t_h_dense),intent(inout)    :: this
    class(t_H_base),intent(in)             :: H_in

    select type(H_in)
    class is(t_h_dense)
        this%H=this%H+H_in%H
    class default
        STOP "Cannot add t_h_dense type with Hamiltonian that is not a class of t_h_dense"
    end select

end subroutine 

subroutine destroy_child(this)
    class(t_h_dense),intent(inout)    :: this

    if(this%is_set())then
        deallocate(this%H)
    endif
end subroutine

subroutine set_from_Hcoo(this,H_coo)
    class(t_h_dense),intent(inout)  :: this
    type(t_H_coo),intent(inout)     :: H_coo

    !local
    integer                 :: nnz,i
    real(8),allocatable     :: val(:)
    integer,allocatable     :: rowind(:),colind(:)

    Call H_coo%pop_par(this%dimH,nnz,val,rowind,colind)
    allocate(this%H(this%dimH(1),this%dimH(2)),source=0.0d0)
    do i=1,nnz
        this%H(rowind(i),colind(i))=val(i)
    enddo
end subroutine 


subroutine eval_single(this,E,i_m,dim_bnd,lat)
    ! input
    class(t_h_dense),intent(in)     :: this
    type(lattice), intent(in)       :: lat
    integer, intent(in)             :: i_m
    integer, intent(in)             :: dim_bnd(2,number_different_order_parameters)    !starting/final index in respective dim_mode of the order parameter (so that energy of single magnetic atom can be be calculated
    ! output
    real(8), intent(out)            :: E
    ! internal
    real(8),pointer                 :: modes_l(:),modes_r(:)
    real(8),allocatable,target      :: vec_l(:),vec_r(:)
    real(8)                         :: tmp(this%dimH(2))
    integer                         :: bnd(2)

    ERROR STOP "THIS PROBABLY NO LONGER WORKS WITH THE NEW MODE_L/MODE_R"   !and in general might be much more difficult to implement with eg. rank 4 in M-space only
    Call lat%point_order(this%op_l,this%dimH(1),modes_l,vec_l)
    Call lat%point_order_single(this%op_r,i_m,dim_bnd,this%dim_mode(2),modes_r,vec_r,bnd)

    tmp=matmul(this%H(:,bnd(1):bnd(2)),modes_r)
    E=dot_product(modes_l,tmp)

    if(allocated(vec_l)) deallocate(vec_l)
    if(allocated(vec_r)) deallocate(vec_r)
end subroutine 


end module
