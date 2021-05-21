module m_H_coo
!Hamiltonian type only for calculating coo parameters without external library
!hence evaluation does not work <- don't use in ham_type_gen 
use m_H_deriv, only: t_H
use m_H_type, only: t_H_base
use m_derived_types, only: lattice, number_different_order_parameters,op_abbrev_to_int
use mpi_basic                

type,extends(t_H) :: t_H_coo
    private
    integer             :: nnz=0  !number of entries in sparse matrix
    real(8),allocatable :: val(:)
    integer,allocatable :: rowind(:),colind(:)
contains
    !necessary t_H routines
    procedure :: eval_single
    procedure :: init_connect
    procedure :: init_mult_connect_2
    procedure :: init_coo

    procedure :: destroy_child
    procedure :: copy_child
    procedure :: add_child

    procedure :: optimize
    procedure :: mult_l,mult_r
    procedure :: mult_l_cont,mult_r_cont
    procedure :: mult_l_disc,mult_r_disc

    !MPI
    procedure :: send
    procedure :: recv
    procedure :: distribute
    procedure :: bcast

    !routine to get all coo parameters 
    !WARNING, DESTROYS INSTANCE
    procedure :: pop_par
end type
private
public t_H,t_H_coo
contains 


subroutine mult_r_cont(this,bnd,vec,res)
    class(t_H_coo),intent(in)    :: this
    integer,intent(in)           :: bnd(2)
    real(8),intent(in)           :: vec(bnd(2)-bnd(1)+1)
    real(8),intent(inout)        :: res(:)   !result matrix-vector product

    STOP "IMPLEMENT if necessary"
end subroutine 

subroutine mult_l_cont(this,bnd,vec,res)
    class(t_H_coo),intent(in)    :: this
    integer,intent(in)           :: bnd(2)
    real(8),intent(in)           :: vec(bnd(2)-bnd(1)+1)
    real(8),intent(inout)        :: res(:)   !result matrix-vector product

    STOP "IMPLEMENT if necessary"
end subroutine 

subroutine mult_r_disc(this,N,ind,vec,res)
    class(t_H_coo),intent(in)    :: this
    integer,intent(in)           :: N
    integer,intent(in)           :: ind(N)
    real(8),intent(in)           :: vec(N)
    real(8),intent(inout)        :: res(:)   !result matrix-vector product

    STOP "IMPLEMENT if necessary"
end subroutine 

subroutine mult_l_disc(this,N,ind,vec,res)
    class(t_H_coo),intent(in)    :: this
    integer,intent(in)           :: N
    integer,intent(in)           :: ind(N)
    real(8),intent(in)           :: vec(N)
    real(8),intent(inout)        :: res(:)   !result matrix-vector product

    STOP "IMPLEMENT if necessary"
end subroutine 


subroutine mult_r(this,lat,res)
    class(t_H_coo),intent(in)    :: this
    type(lattice),intent(in)     :: lat
    real(8),intent(inout)        :: res(:)

    STOP "IMPLEMENT mult_r FOR t_H_coo in m_H_coo if really necessary"
end subroutine 

subroutine mult_l(this,lat,res)
    class(t_H_coo),intent(in)    :: this
    type(lattice),intent(in)     :: lat
    real(8),intent(inout)        :: res(:)

    STOP "IMPLEMENT mult_l FOR t_H_coo in m_H_coo if really necessary"
end subroutine 

subroutine optimize(this)
    class(t_H_coo),intent(inout)    :: this

    STOP "IMPLEMENT optimize FOR t_H_coo in m_H_coo if really necessary"
end subroutine 


subroutine destroy_child(this)
    class(t_H_coo),intent(inout)    :: this

    this%dimH=0
    this%nnz=0
    if(allocated(this%val   )) deallocate(this%val   )
    if(allocated(this%rowind)) deallocate(this%rowind)
    if(allocated(this%colind)) deallocate(this%colind)
end subroutine 

subroutine copy_child(this,Hout)
    class(t_H_coo),intent(in)       :: this
    class(t_H_base),intent(inout)   :: Hout

    select type(Hout)
    class is(t_H_coo)
        Hout%nnz=this%nnz
        allocate(Hout%val   ,source=this%val   )
        allocate(Hout%rowind,source=this%rowind)
        allocate(Hout%colind,source=this%colind)
    class default
        STOP "Cannot copy t_H_coo type with Hamiltonian that is not a class of t_H_coo"
    end select
end subroutine 

subroutine add_child(this,H_in)
    class(t_H_coo),intent(inout)    :: this
    class(t_H_base),intent(in)           :: H_in

    STOP "IMPLEMENT ADDIND FOR t_H_coo in m_H_coo if really necessary"
end subroutine 

subroutine pop_par(this,dimH,nnz,val,rowind,colind)
    class(t_H_coo),intent(inout)    :: this
    integer,intent(out)                 :: dimH(2)
    integer,intent(out)                 :: nnz
    real(8),allocatable,intent(out)     :: val(:)
    integer,allocatable,intent(out)     :: rowind(:),colind(:)

    dimH=this%dimH
    nnz=this%nnz
    this%dimH=0
    this%nnz=0
    Call MOVE_ALLOC(this%val,val)
    Call MOVE_ALLOC(this%rowind,rowind)
    Call MOVE_ALLOC(this%colind,colind)

end subroutine

subroutine init_coo(this,rowind,colind,val,dim_mode,op_l,op_r,lat,mult_M_single)
    !constructs a Hamiltonian based directly on the coo-arrays, which are moved into the type
    class(t_H_coo),intent(inout)        :: this
    real(8),allocatable,intent(inout)   :: val(:)
    integer,allocatable,intent(inout)   :: rowind(:),colind(:)
    integer,intent(in)                  :: dim_mode(2)
    character(len=*),intent(in)         :: op_l         !which order parameters are used at left  side of local Hamiltonian-matrix
    character(len=*),intent(in)         :: op_r         !which order parameters are used at right side of local Hamiltonian-matrix
    type(lattice),intent(in)            :: lat
    integer,intent(in)                  :: mult_M_single !gives the multiple with which the energy_single calculation has to be multiplied (1 for on-site terms, 2 for eg. magnetic exchange)


    integer,allocatable :: order_l(:),order_r(:)

    if(this%is_set()) STOP "cannot set hamiltonian as it is already set"

    allocate(order_l(len(op_l)),source=0)
    allocate(order_r(len(op_r)),source=0)
    order_l=op_abbrev_to_int(op_l)
    order_r=op_abbrev_to_int(op_r)

    if(.not.allocated(rowind)) ERROR STOP "valind has to be allocated to initialize t_H_coo though init_coo"
    if(.not.allocated(colind)) ERROR STOP "rowind has to be allocated to initialize t_H_coo though init_coo"
    if(.not.allocated(val   )) ERROR STOP "val    has to be allocated to initialize t_H_coo though init_coo"

    this%nnz=size(val)

    Call move_alloc(rowind,this%rowind)
    Call move_alloc(colind,this%colind)
    Call move_alloc(val,   this%val   )

    Call this%init_base(lat,order_l,order_r)
    this%dimH=lat%Ncell*dim_mode
    this%mult_M_single=mult_M_single
    Call check_H(this)
end subroutine 


subroutine init_connect(this,connect,Hval,Hval_ind,order,lat,mult_M_single)
    !constructs a Hamiltonian based on only one kind of Hamiltonian subsection (one Hval set)
    class(t_H_coo),intent(inout)    :: this

    type(lattice),intent(in)        :: lat
    character(2),intent(in)         :: order
    real(8),intent(in)              :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
    integer,intent(in)              :: Hval_ind(:,:)
    integer,intent(in)              :: connect(:,:)  !(2,Nentries) index in (1:Ncell) basis of both connected sites 
    integer,intent(in)              :: mult_M_single !gives the multiple with which the energy_single calculation has to be multiplied (1 for on-site terms, 2 for eg. magnetic exchange)

    integer             :: order_int(2)
    integer             :: dim_mode(2)
    integer             :: nnz
    integer             :: i
    integer             :: N_connect,sizeHin

    if(this%is_set()) STOP "cannot set hamiltonian as it is already set"
    N_connect=size(connect,2)
    sizeHin=size(Hval)
    nnz=sizeHin*N_connect

    !fill temporary coordinate format spare matrix
    allocate(this%val(nnz),source=0.0d0)
    allocate(this%colind(nnz),source=0)
    allocate(this%rowind(nnz),source=0)
    order_int=op_abbrev_to_int(order)
    dim_mode(1)=lat%get_order_dim(order_int(1))
    dim_mode(2)=lat%get_order_dim(order_int(2))

    this%nnz=sizeHin*N_connect

!$omp parallel do 
    do i=1,N_connect
        this%rowind(1+(i-1)*sizeHin:i*sizeHin)=(connect(1,i)-1)*dim_mode(1)+Hval_ind(1,:)
        this%colind(1+(i-1)*sizeHin:i*sizeHin)=(connect(2,i)-1)*dim_mode(2)+Hval_ind(2,:)
        this%val   (1+(i-1)*sizeHin:i*sizeHin)=Hval(:)
    enddo
!$omp end parallel do 

    Call this%init_base(lat,[order_int(1)],[order_int(2)])
    this%dimH=lat%Ncell*dim_mode
    this%mult_M_single=mult_M_single
    Call check_H(this)
end subroutine 

subroutine init_mult_connect_2(this,connect,Hval,Hval_ind,op_l,op_r,lat,mult_M_single,dim_mode_in)
    !Constructs a Hamiltonian that depends on more than 2 order parameters but only at 2 sites (i.e. some terms are onsite)
    !(example: ME-coupling M_i*E_i*M_j
    use m_derived_types, only: lattice,op_abbrev_to_int
    class(t_H_coo),intent(inout)    :: this

    type(lattice),intent(in)        :: lat
    !input Hamiltonian
    real(8),intent(in)              :: Hval(:)  !values of local Hamiltonian for each line
    integer,intent(in)              :: Hval_ind(:,:)  !indices in order-parameter space for Hval
    character(len=*),intent(in)     :: op_l         !which order parameters are used at left  side of local Hamiltonian-matrix
    character(len=*),intent(in)     :: op_r         !which order parameters are used at right side of local Hamiltonian-matrix
    integer,intent(in)              :: connect(:,:) !lattice sites to be connected (2,Nconnections)
    integer,intent(in)              :: mult_M_single !gives the multiple with which the energy_single calculation has to be multiplied (1 for on-site terms, 2 for eg. magnetic exchange)
    integer,intent(in),optional     :: dim_mode_in(2)   !optional way of putting in dim_mode directly (mainly for custom(not fully unfolded)rankN tensors)

    integer,allocatable :: order_l(:),order_r(:)
    integer             :: dim_mode(2)
    integer             :: nnz
    integer             :: i,j
    integer             :: ival
    integer             :: N_connect

    if(this%is_set()) STOP "cannot set hamiltonian as it is already set"
    N_connect=size(connect,2)
    nnz=size(Hval)*N_connect

    allocate(order_l(len(op_l)),source=0)
    allocate(order_r(len(op_r)),source=0)
    order_l=op_abbrev_to_int(op_l)
    order_r=op_abbrev_to_int(op_r)

    this%mult_M_single=mult_M_single

    !fill temporary coordinate format spare matrix
    if(present(dim_mode_in))then
        dim_mode=dim_mode_in
    else
        dim_mode=1
        do i=1,size(order_l)
            dim_mode(1)=dim_mode(1)*lat%get_order_dim(order_l(i))
        enddo
        do i=1,size(order_r)
            dim_mode(2)=dim_mode(2)*lat%get_order_dim(order_r(i))
        enddo
    endif

    !set local H
    allocate(this%val(nnz),source=0.0d0)
    allocate(this%colind(nnz),source=0)
    allocate(this%rowind(nnz),source=0)
    nnz=0
    do j=1,N_connect
        do ival=1,size(Hval)
            nnz=nnz+1
            this%rowind(nnz)=(connect(1,j)-1)*dim_mode(1)+Hval_ind(1,ival)
            this%colind(nnz)=(connect(2,j)-1)*dim_mode(2)+Hval_ind(2,ival)
            this%val(nnz)=Hval(ival)
        enddo
    enddo

    !fill additional type parameters
    this%nnz=nnz
    this%dimH=lat%Ncell*dim_mode
    Call this%init_base(lat,order_l,order_r)
    Call check_H(this)
end subroutine 

subroutine check_H(this)
    !some sanity test for this Hamiltonian
    class(t_H_coo),intent(inout)    :: this

    if(maxval(this%rowind)>this%dimH(1)) STOP "H_coo rowind entry larger than dimH"
    if(maxval(this%colind)>this%dimH(2)) STOP "H_coo colind entry larger than dimH"
    if(any(this%rowind<1)) STOP "H_coo rowind entry smaller than 1"
    if(any(this%colind<1)) STOP "H_coo colind entry smaller than 1"

end subroutine

subroutine eval_single(this,E,i_m,order,lat,work)
    use m_work_ham_single
    ! input
    class(t_H_coo),intent(in)   :: this
    type(lattice), intent(in)   :: lat
    integer, intent(in)         :: i_m
    integer,intent(in)          :: order
    ! output
    real(kind=8), intent(out)    :: E
    type(work_ham_single),intent(inout) :: work

    ERROR STOP "CANNOT calculate single energy contributions with this Hamiltonian-type"
end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!            MPI ROUTINES           !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine send(this,ithread,tag,com)
    class(t_H_coo),intent(in)      :: this
    integer,intent(in)              :: ithread
    integer,intent(in)              :: tag
    integer,intent(in)              :: com

#ifdef CPP_MPI
    Call this%send_base(ithread,tag,com)
    ERROR STOP "IMPLEMENT if really necessary"
#else
    continue
#endif
end subroutine

subroutine recv(this,ithread,tag,com)
    class(t_H_coo),intent(inout)   :: this
    integer,intent(in)              :: ithread
    integer,intent(in)              :: tag
    integer,intent(in)              :: com

#ifdef CPP_MPI
    Call this%recv_base(ithread,tag,com)
    ERROR STOP "IMPLEMENT if really necessary"
#else
    continue
#endif
end subroutine

subroutine bcast(this,comm)
    class(t_H_coo),intent(inout)        ::  this
    type(mpi_type),intent(in)       ::  comm

    STOP "IMPLEMENT bcast FOR t_H_coo in m_H_coo if really necessary"
end subroutine 

subroutine distribute(this,comm)
    class(t_H_coo),intent(inout)        ::  this
    type(mpi_type),intent(in)       ::  comm

    STOP "IMPLEMENT disbribute FOR t_H_coo in m_H_coo if really necessary"
end subroutine 


end module
