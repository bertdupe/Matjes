module m_H_tb_coo
use m_H_tb_base
use m_ham_init_type, only: parameters_ham_init 
implicit none
private
public H_TB_coo, H_TB_coo_based

type,abstract,extends(H_TB) :: H_TB_coo_based
    !type used to avoid having to manually implement set_from_Hcoo for each implementation that is constructed using H_tb_coo
contains
    procedure(int_set_from_Hcoo),deferred   :: set_from_Hcoo    !routine to set the local parameters of the Hamiltonian based from a set H_tb_coo-type
    procedure :: init_connect => init_connect_based 
    procedure :: init_coo     => init_coo_based 
end type

type,extends(H_TB) :: H_tb_coo
    !type implementing easy access Hamiltonian description in coo-sparse format
    !used for initialization of more specified Hamiltonian desciptions
    private
    integer                 :: nnz=0  !number of entries in sparse matrix
    complex(8),allocatable  :: val(:)
    integer,allocatable     :: rowind(:),colind(:)
contains
    procedure   :: init_coo
    procedure   :: init_connect
    procedure   :: add_child
    procedure   :: copy_child
    procedure   :: destroy_child
    procedure   :: pop_par
    procedure   :: get_par
    procedure   :: mv

    procedure   :: to_BdG

    procedure   :: mult_r
    procedure   :: get_eval
    procedure   :: get_evec
end type

abstract interface
    subroutine int_set_from_Hcoo(this,H_coo)
        import H_tb_coo,H_TB_coo_based
        class(H_TB_coo_based),intent(inout)  :: this
        type(H_tb_coo),intent(inout)       :: H_coo
    end subroutine
end interface


contains

subroutine to_BdG(this)
    !adds the BdG hole/hole part if only the electronic part is described
    class(H_TB_coo),intent(inout)   :: this
    integer                     :: nnz
    complex(8),allocatable      :: val(:)
    integer,allocatable         :: rowind(:),colind(:)
    type(parameters_ham_init)   :: hinit
    integer                     :: nBdG,i

    if(this%nsc==1)then
        Call this%get_hinit(hinit)
        nBdG=hinit%norb*hinit%nspin*hinit%ncell
        hinit%nsc=2
        Call this%pop_par(nnz,val,rowind,colind)
        Call this%init_base(hinit)
        this%nnz=2*nnz
        allocate(this%val(this%nnz),this%rowind(this%nnz),this%colind(this%nnz))
        this%val   (1:nnz)=val  
        this%rowind(1:nnz)=rowind
        this%colind(1:nnz)=colind
        this%val   (nnz+1:2*nnz)=conjg(val)
        this%rowind(nnz+1:2*nnz)=rowind+nBdG
        this%colind(nnz+1:2*nnz)=colind+nBdG
        forall(i=1+nnz:2*nnz,modulo(this%rowind(i),2)==modulo(this%colind(i),2)) this%val(i)=-this%val(i)   !minus for spin-conserving terms
    elseif(this%nsc/=2)then
        STOP "Unexected nsc"
    endif

end subroutine


subroutine mv(this,Hout)
    class(H_TB_coo),intent(inout)   :: this
    class(H_TB),intent(inout)       :: Hout
    
    select type(Hout)
    type is(H_TB_coo)
        Call Hout%init_otherH(this)
        Hout%nnz=this%nnz
        Call move_alloc(this%val,Hout%val)
        Call move_alloc(this%rowind,Hout%rowind)
        Call move_alloc(this%colind,Hout%colind)
    class is(H_TB_coo_based)
        Call Hout%set_from_Hcoo(this)
    class default
        STOP "Cannot move H_TB_coo type with Hamiltonian that is not a class of H_TB_coo or H_TB_coo_based"
    end select
    Call this%destroy()
end subroutine


subroutine get_par(this,val,rowind,colind)
    class(H_TB_coo),intent(in)          :: this
    complex(8),allocatable,intent(out)  :: val(:)
    integer,allocatable,intent(out)     :: rowind(:),colind(:)

    val=this%val
    rowind=this%rowind
    colind=this%colind
end subroutine

subroutine pop_par(this,nnz,val,rowind,colind)
    class(H_TB_coo),intent(inout)     :: this
    integer,intent(out)                 :: nnz
    complex(8),allocatable,intent(out)  :: val(:)
    integer,allocatable,intent(out)     :: rowind(:),colind(:)

    nnz=this%nnz
    this%dimH=0
    this%nnz=0
    Call MOVE_ALLOC(this%val,val)
    Call MOVE_ALLOC(this%rowind,rowind)
    Call MOVE_ALLOC(this%colind,colind)
    Call this%destroy
end subroutine


subroutine init_connect_based(this,connect,Hval,Hval_ind,io,ind_offset)
    use m_derived_types, only: lattice
    class(H_TB_coo_based),intent(inout)     :: this
    complex(8),intent(in)                   :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
    integer,intent(in)                      :: Hval_ind(:,:)
    integer,intent(in)                      :: connect(:,:)
    type(parameters_ham_init),intent(in)    :: io
    integer,intent(in),optional             :: ind_offset(2)
    !local
    type(H_TB_coo)    :: H_coo

    Call H_coo%init_connect(connect,Hval,Hval_ind,io,ind_offset)
    Call this%set_from_Hcoo(H_coo)
end subroutine 


subroutine init_coo_based(this,val,row,col,io,ind_offset)
    use m_derived_types, only: lattice
    class(H_TB_coo_based),intent(inout)     :: this
    complex(8),intent(in)                   :: val(:)  !all entries between 2 cell sites of considered orderparameter
    integer,intent(in)                      :: row(size(val)),col(size(val))
    type(parameters_ham_init),intent(in)    :: io
    integer,intent(in),optional             :: ind_offset(2)
    !local
    type(H_TB_coo)    :: H_coo

    Call H_coo%init_coo(val,row,col,io,ind_offset)
    Call this%set_from_Hcoo(H_coo)
end subroutine 


subroutine init_coo(this,val,row,col,io,ind_offset)
    !constructs Hamiltonian manually for a given sparse coo-input
    class(H_tb_coo),intent(inout)           :: this
    complex(8),intent(in)                   :: val(:)  !all entries between 2 cell sites of considered orderparameter
    integer,intent(in)                      :: row(size(val)),col(size(val))
    type(parameters_ham_init),intent(in)    :: io
    integer,intent(in),optional             :: ind_offset(2) !integer offset in each Hamiltonian dimension (to easily access other sectors in BdG-space)

    if(this%is_set()) STOP "cannot set hamiltonian as it is already set"
    Call this%init_base(io)
    this%nnz=size(val)
    this%val=val
    this%rowind=row
    this%colind=col
    if(present(ind_offset))then
        this%rowind=this%rowind+ind_offset(1)
        this%colind=this%colind+ind_offset(2)
    endif
end subroutine


subroutine init_connect(this,connect,Hval,Hval_ind,io,ind_offset)
    !constructs a Hamiltonian based on only one kind of Hamiltonian subsection (one Hval set)
    !takes the Hamiltonian entries within one unit-cell (Hval, H_ind) and unfolds them on the supercell using the connections
    !for one difference vector given by the input: connect
    class(H_tb_coo),intent(inout)           :: this             !Hamiltonian type to be set
    complex(8),intent(in)                   :: Hval(:)          !all entries between 2 cell sites of considered orderparameter
    integer,intent(in)                      :: Hval_ind(:,:)    !(2,sizeHin) indices (column, row) within unit-cell
    integer,intent(in)                      :: connect(:,:)     !(2,Nentries) index in (1:Ncell) basis of both connected sites (super-cell)
    type(parameters_ham_init),intent(in)    :: io               !shape/dimension parameters to be set
    integer,intent(in),optional             :: ind_offset(2)    !integer offset in each Hamiltonian dimension (to easily access other sectors in BdG-space)

    integer :: nnz
    integer :: i
    integer :: N_connect,sizeHin

    if(this%is_set()) STOP "cannot set hamiltonian as it is already set"
    Call this%init_base(io) !initialize correct size
    N_connect=size(connect,2)
    sizeHin=size(Hval)
    nnz=sizeHin*N_connect

    !fill temporary coordinate format spare matrix
    allocate(this%val(nnz),source=(0.0d0,0.0d0))
    allocate(this%colind(nnz),source=0)
    allocate(this%rowind(nnz),source=0)

    this%nnz=sizeHin*N_connect

    if(present(ind_offset))then
!$omp parallel do 
        do i=1,N_connect
            this%rowind(1+(i-1)*sizeHin:i*sizeHin)=(connect(1,i)-1)*this%ndim+Hval_ind(1,:)+ind_offset(1)
            this%colind(1+(i-1)*sizeHin:i*sizeHin)=(connect(2,i)-1)*this%ndim+Hval_ind(2,:)+ind_offset(2)
            this%val   (1+(i-1)*sizeHin:i*sizeHin)=Hval(:)
        enddo
!$omp end parallel do 
    else
!$omp parallel do 
        do i=1,N_connect
            this%rowind(1+(i-1)*sizeHin:i*sizeHin)=(connect(1,i)-1)*this%ndim+Hval_ind(1,:)
            this%colind(1+(i-1)*sizeHin:i*sizeHin)=(connect(2,i)-1)*this%ndim+Hval_ind(2,:)
            this%val   (1+(i-1)*sizeHin:i*sizeHin)=Hval(:)
        enddo
!$omp end parallel do 
    endif
end subroutine 

subroutine add_child(this,H_in)
    class(H_tb_coo),intent(inout)     :: this
    class(H_tb),intent(in)              :: H_in
    !tmp
    complex(8),allocatable :: val(:)
    integer,allocatable :: rowind(:),colind(:)

    select type(H_in)
    class is(H_tb_coo)
        Call move_alloc(this%val   ,val   )
        Call move_alloc(this%rowind,rowind)
        Call move_alloc(this%colind,colind)
        this%nnz=this%nnz+H_in%nnz
        this%val   =[val   ,H_in%val   ]
        this%rowind=[rowind,H_in%rowind]
        this%colind=[colind,H_in%colind]
    class default
        STOP "Cannot add H_TB_coo type with Hamiltonian that is not a class of H_TB_coo"
    end select
end subroutine

recursive subroutine copy_child(this,Hout)
    class(H_TB_coo),intent(in)  :: this
    class(H_TB),intent(inout)   :: Hout
    type(H_TB_coo)              :: Htmp
    
    select type(Hout)
    type is(H_TB_coo)
        Hout%nnz=this%nnz
        allocate(Hout%val,source=this%val)
        allocate(Hout%rowind,source=this%rowind)
        allocate(Hout%colind,source=this%colind)
    class is(H_TB_coo_based)
        Call this%copy(Htmp)
        Call Hout%set_from_Hcoo(Htmp)
        Call Htmp%destroy()
    class default
        STOP "Cannot copy H_TB_coo type with Hamiltonian that is not a class of H_TB_coo"
    end select
end subroutine

subroutine destroy_child(this)
    class(H_tb_coo),intent(inout)    :: this

    if(this%is_set())then
        this%nnz=0
        if(allocated(this%val   )) deallocate(this%val   )
        if(allocated(this%colind)) deallocate(this%colind)
        if(allocated(this%rowind)) deallocate(this%rowind)
    endif
end subroutine


subroutine mult_r(this,vec,res,alpha,beta)
    class(H_tb_coo),intent(in)      :: this
    complex(8),intent(in   )        :: vec(this%dimH)
    complex(8),intent(inout)        :: res(this%dimH)
    complex(8),intent(in),optional     :: alpha
    complex(8),intent(in),optional     :: beta
    ERROR STOP "Multiply_r from H_TB_coo type is not implemented"
end subroutine

subroutine get_eval(this,eval)
    class(H_TB_coo),intent(in)      ::  this
    real(8),intent(out),allocatable ::  eval(:)

    ERROR STOP "Get eigenvalues from H_TB_coo type is not implemented"
    allocate(eval(1),source=0.0d0)
end subroutine

subroutine get_evec(this,eval,evec)
    class(H_TB_coo),intent(in)  ::  this
    real(8),allocatable,intent(out)         ::  eval(:)
    complex(8),allocatable,intent(out)      ::  evec(:,:)

    ERROR STOP "Get eigenvectors from H_TB_coo type is not implemented"
    allocate(eval(1),source=0.0d0)
    allocate(evec(1,1),source=(0.0d0,0.0d0))
end subroutine


end module
