module m_H_type_coo
!Hamiltonian type only for calculating coo parameters without external library
!hence evaluation does not work <- don't use in ham_type_gen 
use m_H_type

type,extends(t_H) :: t_H_coo
    private
    integer             :: dimH(2)=0 !dimension of Hamiltonian
    integer             :: nnz=0  !number of entries in sparse matrix
    real(8),allocatable :: val(:)
    integer,allocatable :: rowind(:),colind(:)
contains
    !necessary t_H routines
    procedure :: eval_single
    procedure :: eval_all
    procedure :: set_H
    procedure :: set_H_1
    procedure :: add_H
    procedure :: destroy
    procedure :: optimize
    procedure :: copy
    procedure :: mult_l,mult_r

    !routine to get all coo parameters 
    !WARNING, DESTROYS INSTANCE
    procedure :: pop_par
end type
private
public t_H,t_H_coo
contains 

subroutine mult_r(this,lat,vec)
    use m_derived_types, only: lattice
    class(t_H_coo),intent(in)    :: this
    type(lattice),intent(in)     :: lat
    real(8),intent(inout)        :: vec(:)

    STOP "IMPLEMENT apply_r FOR t_H_coo in m_H_type_coo if really necessary"
end subroutine 

subroutine mult_l(this,lat,vec)
    use m_derived_types, only: lattice
    class(t_H_coo),intent(in)    :: this
    type(lattice),intent(in)     :: lat
    real(8),intent(inout)        :: vec(:)

    STOP "IMPLEMENT apply_l FOR t_H_coo in m_H_type_coo if really necessary"
end subroutine 


subroutine optimize(this)
    class(t_H_coo),intent(inout)    :: this

    STOP "IMPLEMENT optimize FOR t_H_coo in m_H_type_coo if really necessary"
end subroutine 

subroutine destroy(this)
    class(t_H_coo),intent(inout)    :: this

    !inherited data
    Call this%set_prepared(.false.)
    deallocate(this%op_l,this%op_r)

    this%dimH=0
    this%nnz=0
    deallocate(this%val,this%rowind,this%colind)
end subroutine 

subroutine copy(this,Hout)
    class(t_H_coo),intent(in)    :: this
    class(t_H),intent(inout)     :: Hout

    STOP "IMPLEMENT copy FOR t_H_coo in m_H_type_coo if really necessary"

end subroutine 

subroutine add_H(this,H_add)
    class(t_H_coo),intent(inout)    :: this
    class(t_H),intent(in)           :: H_add

    STOP "IMPLEMENT ADDIND FOR t_H_coo in m_H_type_coo if really necessary"

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


subroutine set_H_1(this,line,Hval,Hval_ind,order,lat)
    !constructs a Hamiltonian based on only one kind of Hamiltonian subsection (one Hval set)
    use m_derived_types, only: operator_real_order_N,lattice
    class(t_H_coo),intent(inout)    :: this

    type(lattice),intent(in)        :: lat
    integer,intent(in)              :: order(2)
    real(8),intent(in)              :: Hval(:)  !all entries between 2 cell sites of considered orderparameter
    integer,intent(in)              :: Hval_ind(:,:)
    integer,intent(in)              :: line(:,:)

    integer             :: dim_mode(2)
    integer             :: nnz
    integer             :: i,j,l,ii
    integer             :: i1,i2,ival
    integer             :: N_neigh
    integer             :: N_persite,N_site

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

    ii=0
    do i=1,N_site
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
    allocate(this%op_l(1),source=order(1))
    allocate(this%op_r(1),source=order(2))
    this%nnz=ii
    this%dimH=N_site*dim_mode
    allocate(this%colind,source=colind(1:this%nnz))
    deallocate(colind)
    allocate(this%rowind,source=rowind(1:this%nnz))
    deallocate(rowind)
    allocate(this%val,source=val(1:this%nnz))
    deallocate(val)
    Call this%set_prepared(.true.)

end subroutine 


subroutine set_H(this,energy_in,lat)
    use m_derived_types, only: operator_real_order_N,lattice
    class(t_H_coo),intent(inout)  :: this
    type(operator_real_order_N)     :: energy_in
    type(lattice),intent(in)        :: lat

    integer             :: nnz
    integer             :: i,j,l,ii
    integer             :: i1,i2
    integer             :: N_neigh
    integer             :: N_persite,N_site

    real(8),allocatable :: val(:)
    integer,allocatable :: rowind(:),colind(:)


    N_neigh=size(energy_in%line,1)
    N_site=size(energy_in%line,2)
    
    !estimate size
    N_persite=0
    do i=1,N_neigh
        N_persite=N_persite+count(energy_in%value(i,1)%order_op(1)%Op_loc /= 0.0d0)
    enddo
    nnz=N_persite*N_site

    !fill temporary coordinate format spare matrix
    allocate(val(nnz),source=0.0d0)
    allocate(colind(nnz),source=0)
    allocate(rowind(nnz),source=0)
    ii=0
    do i=1,N_site
        do l=1,N_neigh
            j=energy_in%line(l,i)
            do i2=1,lat%dim_mode
                do i1=1,lat%dim_mode
                    if(energy_in%value(l,i)%order_op(1)%Op_loc(i1,i2)/= 0.0d0)then
                        ii=ii+1
                        colind(ii)=(j-1)*lat%dim_mode+i1
                        rowind(ii)=(i-1)*lat%dim_mode+i2
                        val(ii)=energy_in%value(l,i)%order_op(1)%Op_loc(i1,i2)
                    endif
                enddo
            enddo
        enddo
    enddo

    !fill type
    this%nnz=ii
    this%dimH=N_site*lat%dim_mode
    allocate(this%colind,source=colind(1:this%nnz))
    deallocate(colind)
    allocate(this%rowind,source=rowind(1:this%nnz))
    deallocate(rowind)
    allocate(this%val,source=val(1:this%nnz))
    deallocate(val)
    Call this%set_prepared(.true.)

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


subroutine eval_all(this,E,lat)
    use m_derived_types, only: lattice
    class(t_H_coo),intent(in)    :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(out)            :: E

    STOP "CANNOT EVALUATE t_H_coo"
    !alternatively add some evaluation without an library
    
end subroutine 

end module
