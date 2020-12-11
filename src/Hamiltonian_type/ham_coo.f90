module m_H_coo
!Hamiltonian type only for calculating coo parameters without external library
!hence evaluation does not work <- don't use in ham_type_gen 
use m_H_type

type,extends(t_H) :: t_H_coo
    private
    integer             :: nnz=0  !number of entries in sparse matrix
    real(8),allocatable :: val(:)
    integer,allocatable :: rowind(:),colind(:)
contains
    !necessary t_H routines
    procedure :: eval_single
    procedure :: init_1
    procedure :: init_mult_2

    procedure :: destroy_child
    procedure :: copy_child
    procedure :: add_child

    procedure :: optimize
    procedure :: mult_l,mult_r

    !routine to get all coo parameters 
    !WARNING, DESTROYS INSTANCE
    procedure :: pop_par
end type
private
public t_H,t_H_coo
contains 

subroutine mult_r(this,lat,res)
    use m_derived_types, only: lattice
    class(t_H_coo),intent(in)    :: this
    type(lattice),intent(in)     :: lat
    real(8),intent(inout)        :: res(:)

    STOP "IMPLEMENT mult_r FOR t_H_coo in m_H_coo if really necessary"
end subroutine 

subroutine mult_l(this,lat,res)
    use m_derived_types, only: lattice
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
    deallocate(this%val,this%rowind,this%colind)
end subroutine 

subroutine copy_child(this,Hout)
    class(t_H_coo),intent(in)    :: this
    class(t_H),intent(inout)     :: Hout

    STOP "IMPLEMENT copy FOR t_H_coo in m_H_coo if really necessary"

end subroutine 

subroutine add_child(this,H_in)
    class(t_H_coo),intent(inout)    :: this
    class(t_H),intent(in)           :: H_in

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


subroutine init_1(this,line,Hval,Hval_ind,order,lat,mult_M_single)
    !constructs a Hamiltonian based on only one kind of Hamiltonian subsection (one Hval set)
    use m_derived_types, only: lattice
    class(t_H_coo),intent(inout)    :: this

    type(lattice),intent(in)        :: lat
    integer,intent(in)              :: order(2)
    real(8),intent(in)              :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
    integer,intent(in)              :: Hval_ind(:,:)
    integer,intent(in)              :: line(:,:)
    integer,intent(in)              :: mult_M_single !gives the multiple with which the energy_single calculation has to be multiplied (1 for on-site terms, 2 for eg. magnetic exchange)

    integer             :: dim_mode(2)
    integer             :: nnz
    integer             :: i,j,l,ii
    integer             :: ival
    integer             :: N_neigh,N_site

    real(8),allocatable :: val(:)
    integer,allocatable :: rowind(:),colind(:)


    N_neigh=size(line,1)
    N_site=size(line,2)
    
    nnz=size(Hval)*N_site*N_neigh

    !fill temporary coordinate format spare matrix
    allocate(val(nnz),source=0.0d0)
    allocate(colind(nnz),source=0)
    allocate(rowind(nnz),source=0)
    dim_mode(1)=lat%get_order_dim(order(1))
    dim_mode(2)=lat%get_order_dim(order(2))

    this%mult_M_single=mult_M_single

    ii=0
    do i=1,N_site
        !this could give problems where there are different numbers of neighborsper line
        do l=1,N_neigh
            j=line(l,i)
            do ival=1,size(Hval) 
                ii=ii+1
                colind(ii)=(j-1)*dim_mode(1)+Hval_ind(1,ival)
                rowind(ii)=(i-1)*dim_mode(2)+Hval_ind(2,ival)
                val(ii)=Hval(ival)
            enddo
        enddo
    enddo

    !fill type
    this%nnz=ii
    this%dimH=N_site*dim_mode
    allocate(this%colind,source=colind(1:this%nnz))
    deallocate(colind)
    allocate(this%rowind,source=rowind(1:this%nnz))
    deallocate(rowind)
    allocate(this%val,source=val(1:this%nnz))
    deallocate(val)
    Call this%init_base(lat,order(1:1),order(2:2))
    Call check_H(this)

end subroutine 


subroutine init_mult_2(this,connect,Hval,Hval_ind,op_l,op_r,lat,mult_M_single)
    !Constructs a Hamiltonian that depends on more than 2 order parameters but only at 2 sites (i.e. some terms are onsite)
    !(example: ME-coupling M_i*E_i*M_j
    use m_derived_types, only: lattice
    class(t_H_coo),intent(inout)    :: this

    type(lattice),intent(in)        :: lat
    !input Hamiltonian
    real(8),intent(in)              :: Hval(:)  !values of local Hamiltonian for each line
    integer,intent(in)              :: Hval_ind(:,:)  !indices in order-parameter space for Hval
    integer,intent(in)              :: op_l(:),op_r(:) !which order parameters are used at left/right side of local Hamiltonian-matrix
    integer,intent(in)              :: connect(:,:) !lattice sites to be connected (2,Nconnections)
    integer,intent(in)              :: mult_M_single !gives the multiple with which the energy_single calculation has to be multiplied (1 for on-site terms, 2 for eg. magnetic exchange)

    integer             :: dim_mode(2)
    integer             :: nnz
    integer             :: i,j
    integer             :: ival
    integer             :: N_connect

    N_connect=size(connect,2)
    nnz=size(Hval)*N_connect

    !fill temporary coordinate format spare matrix
    dim_mode=1
    do i=1,size(op_l)
        dim_mode(1)=dim_mode(1)*lat%get_order_dim(op_l(i))
    enddo
    do i=1,size(op_r)
        dim_mode(2)=dim_mode(2)*lat%get_order_dim(op_r(i))
    enddo


    this%mult_M_single=mult_M_single

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
    Call this%init_base(lat,op_l,op_r)
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

subroutine eval_single(this,E,i_m,lat)
    use m_derived_types, only: lattice
    ! input
    class(t_H_coo),intent(in)    :: this
    type(lattice), intent(in)       :: lat
    integer, intent(in)             :: i_m
    ! output
    real(kind=8), intent(out)       :: E

    STOP "CANNOT EVALUATE t_H_coo"
    !alternatively add some evaluation without an library

end subroutine 

end module
