module m_H_dense
!Hamiltonian type specifications using dense matrices and no external library
use m_H_type
use m_H_coo
use m_derived_types, only: lattice


type,extends(t_H) :: t_H_dense
    real(8),allocatable   :: H(:,:)
contains
    !necessary t_H routines
    procedure :: eval_single
    procedure :: init_1    
    procedure :: init_connect    
    procedure :: init_mult_connect_2
    procedure :: init_mult_2   

    procedure :: add_child 
    procedure :: bcast_child 
    procedure :: destroy_child    
    procedure :: copy_child 

    procedure :: optimize
    procedure :: mult_r,mult_l
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

    Call lat%point_order(this%op_r,this%dimH(2),modes,vec)
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

    Call lat%point_order(this%op_l,this%dimH(1),modes,vec)

    if(size(res)/=this%dimH(2)) STOP "size of vec is wrong"
    res=matmul(modes,this%H)
    if(allocated(vec)) deallocate(vec)
    
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

    Call lat%point_order(this%op_r,this%dimH(2),modes,vec)
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

    Call lat%point_order(this%op_l,this%dimH(1),modes,vec)
    if(size(res)/=this%dimH(1)) STOP "size of vec is wrong"
    bnd(1)=this%dim_mode(2)*(i_site-1)+1
    bnd(2)=this%dim_mode(2)*(i_site)
    !terrible implementation striding-wise
    res=matmul(modes,this%H(:,bnd(1):bnd(2)))

    if(allocated(vec)) deallocate(vec)
end subroutine 



subroutine optimize(this)
    class(t_h_dense),intent(inout)   :: this

    !nothing to optimize here
    continue 
end subroutine

subroutine init_connect(this,connect,Hval,Hval_ind,order,lat,mult_M_single)
    use m_derived_types, only: lattice
    class(t_h_dense),intent(inout)    :: this
    type(lattice),intent(in)        :: lat
    character(2),intent(in)         :: order
    real(8),intent(in)              :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
    integer,intent(in)              :: Hval_ind(:,:)
    integer,intent(in)              :: connect(:,:)
    integer,intent(in)              :: mult_M_single

    !local
    type(t_H_coo)           :: H_coo

    Call H_coo%init_connect(connect,Hval,Hval_ind,order,lat,mult_M_single)
    Call set_from_Hcoo(this,H_coo,lat)
end subroutine 


subroutine init_1(this,line,Hval,Hval_ind,order,lat,mult_M_single)
    use m_derived_types, only: lattice
    class(t_h_dense),intent(inout)    :: this

    type(lattice),intent(in)        :: lat
    integer,intent(in)              :: order(2)
    real(8),intent(in)              :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
    integer,intent(in)              :: Hval_ind(:,:)
    integer,intent(in)              :: line(:,:)
    integer,intent(in)              :: mult_M_single

    !local
    type(t_H_coo)           :: H_coo

    Call H_coo%init_1(line,Hval,Hval_ind,order,lat,mult_M_single)
    Call set_from_Hcoo(this,H_coo,lat)
end subroutine 


subroutine init_mult_2(this,connect,Hval,Hval_ind,op_l,op_r,lat,mult_M_single)
    use m_derived_types, only: lattice
    class(t_h_dense),intent(inout)    :: this

    type(lattice),intent(in)        :: lat
    integer,intent(in)              :: op_l(:),op_r(:)
    real(8),intent(in)              :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
    integer,intent(in)              :: Hval_ind(:,:)
    integer,intent(in)              :: connect(:,:)
    integer,intent(in)              :: mult_M_single

    !local
    type(t_H_coo)           :: H_coo

    Call H_coo%init_mult_2(connect,Hval,Hval_ind,op_l,op_r,lat,mult_M_single)
    Call set_from_Hcoo(this,H_coo,lat)
end subroutine 

subroutine init_mult_connect_2(this,connect,Hval,Hval_ind,op_l,op_r,lat,mult_M_single)
    !Constructs a Hamiltonian that depends on more than 2 order parameters but only at 2 sites (i.e. some terms are onsite)
    !(example: ME-coupling M_i*E_i*M_j
    use m_derived_types, only: lattice,op_abbrev_to_int
    class(t_H_dense),intent(inout)    :: this

    type(lattice),intent(in)        :: lat
    !input Hamiltonian
    real(8),intent(in)              :: Hval(:)  !values of local Hamiltonian for each line
    integer,intent(in)              :: Hval_ind(:,:)  !indices in order-parameter space for Hval
    character(len=*),intent(in)     :: op_l         !which order parameters are used at left  side of local Hamiltonian-matrix
    character(len=*),intent(in)     :: op_r         !which order parameters are used at right side of local Hamiltonian-matrix
    integer,intent(in)              :: connect(:,:) !lattice sites to be connected (2,Nconnections)
    integer,intent(in)              :: mult_M_single
    !local
    type(t_H_coo)           :: H_coo

    Call H_coo%init_mult_connect_2(connect,Hval,Hval_ind,op_l,op_r,lat,mult_M_single)
    Call set_from_Hcoo(this,H_coo,lat)
end subroutine

subroutine copy_child(this,Hout)
    class(t_h_dense),intent(in)   :: this
    class(t_H),intent(inout)        :: Hout
    
    select type(Hout)
    class is(t_h_dense)
        allocate(Hout%H,source=this%H)
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
    class(t_H),intent(in)             :: H_in

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

subroutine set_from_Hcoo(this,H_coo,lat)
    type(t_H_coo),intent(inout)         :: H_coo
    type(t_h_dense),intent(inout)     :: this
    type(lattice),intent(in)            :: lat

    !local
    integer                 :: nnz,i
    real(8),allocatable     :: val(:)
    integer,allocatable     :: rowind(:),colind(:)

    Call H_coo%pop_par(this%dimH,nnz,val,rowind,colind)
    allocate(this%H(this%dimH(1),this%dimH(2)),source=0.0d0)
    do i=1,nnz
        this%H(rowind(i),colind(i))=val(i)
    enddo

    Call this%init_base(lat,H_coo%op_l,H_coo%op_r)
    this%mult_M_single=H_coo%mult_M_single
end subroutine 


subroutine eval_single(this,E,i_m,lat)
    use m_derived_types, only: lattice
    ! input
    class(t_h_dense),intent(in)     :: this
    type(lattice), intent(in)       :: lat
    integer, intent(in)             :: i_m
    ! output
    real(kind=8), intent(out)       :: E
    ! internal
    real(8),pointer                 :: modes_l(:),modes_r(:)
    real(8),allocatable,target      :: vec_l(:),vec_r(:)
    real(8)                         :: tmp(this%dimH(2))

    Call lat%point_order(this%op_l,this%dimH(1),modes_l,vec_l)
    Call lat%point_order(this%op_r,this%dimH(2),modes_r,vec_r)

    tmp=matmul(this%H(:,1+(i_m-1)*this%dim_mode(2):i_m*this%dim_mode(2)),modes_r(1+(i_m-1)*this%dim_mode(2):i_m*this%dim_mode(2)))
    E=dot_product(modes_l,tmp)
end subroutine 

end module
