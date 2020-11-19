module m_type_lattice
use m_basic_types
use m_order_parameter
use m_cell,only : t_cell

integer,parameter :: number_different_order_parameters=4    !m,E,B,T
character(len=1),parameter :: order_parameter_abbrev(number_different_order_parameters)=["M","E","B","T"]
character(len=*),parameter :: order_parameter_name(number_different_order_parameters)=[&
                                &   'magnetic   ',&
                                &   'Efield     ',&
                                &   'Bfield     ',&
                                &   'temperature']

! variable that defines the lattice
type lattice
     real(8)        :: areal(3,3) !real space lattice vector
     real(8)        :: astar(3,3) !reciprocal space lattice vector
     type(t_cell)   :: cell
     integer        :: dim_lat(3) !number of unit-cells in each direction in supercell
     integer        :: Ncell !number of unit-cells (product of dim_lat)
     real(8)        :: a_sc(3,3) !real space lattice vectors of supercell
     integer        :: dim_mode !probably obsolete
     integer        :: n_system !no clue what this does
     integer        :: nmag  !this has to be set somewhere consistently
     integer, allocatable :: world(:)
     logical :: periodic(3) !lattice periodic in direction

     real(8),allocatable :: sc_vec_period(:,:)   !real space vectors used to minimize distance using supercell periodicity

!order parameters
     type(order_par)  :: M !magnetization 
     type(order_par)  :: E !Electric field
     type(order_par)  :: B !Magnetic field
     type(order_par)  :: T !Temperature
     !ultimately delete ordpar
     type(order_par)   :: ordpar !obsolete

contains
    !initialization function
    procedure :: init_geo   !initialize the geometric properties (areal,Ncell,...)
    procedure :: init_order !initialize the order parameters dimensions
    !basic function
    procedure :: copy => copy_lattice
    procedure :: copy_val_to => copy_val_lattice
    !get correct order parameter (or combination thereof)
    procedure :: set_order_point
    procedure :: set_order_comb
    procedure :: set_order_comb_exc
    procedure :: set_order_comb_exc_single
    procedure :: point_order => point_order_onsite
    procedure :: point_order_single => point_order_onsite_single
    !!reduce that order parameter again
    procedure :: reduce
    procedure :: reduce_single
    !!reduce vector in order parameter space to the site(e.g. energy Ncell resolution )
    procedure :: reduce_cell
    !real space position functions
    procedure :: get_pos_mag => lattice_get_position_mag
    procedure :: pos_diff_ind => lattice_position_difference
    procedure :: pos_ind => lattice_position
    !index routines
    procedure :: index_m_1 
    procedure :: index_4_1 
    procedure :: index_1_3 
    !misc
    procedure :: used_order => lattice_used_order
    procedure :: get_order_dim 
end type lattice

type point_arr
    !private type used only here
    real(8),pointer     ::  v(:)
end type

private
public lattice, number_different_order_parameters,op_name_to_int,op_int_to_abbrev, t_cell, order_par

contains 

subroutine init_geo(this,areal,alat,dim_lat,boundary)
    !initialize the lattice with the geometric inputs
    use m_vector, only : norm,cross
    use m_constants, only : pi
    class(lattice),intent(inout)    ::  this
    real(8),intent(in)              ::  areal(3,3),alat(3)
    integer,intent(in)              ::  dim_lat(3)
    logical,intent(in)              ::  boundary(3)
    !internal
    real(8)                         :: volume
    integer                         :: i
    !use supercell parameter
    integer,allocatable,dimension(:)     :: i1,i2,i3 
    real(8)         :: tmp_vec(3,3)
    integer         :: j,l,i_vec

    do i=1,3
       this%areal(i,:)=areal(i,:)*alat(i)
    enddo
    this%dim_lat=dim_lat
    this%ncell=product(dim_lat)
    do i=1,3
       this%a_sc(i,:)=this%areal(i,:)*dim_lat(i)
    enddo

    this%periodic=boundary

    ! build up the reciprocal lattice vectors
    volume=dot_product(this%areal(1,:),cross(this%areal(2,:),this%areal(3,:)))
    this%astar(1,:) = cross(this%areal(2,:),this%areal(3,:))/volume
    this%astar(2,:) = cross(this%areal(3,:),this%areal(1,:))/volume
    this%astar(3,:) = cross(this%areal(1,:),this%areal(2,:))/volume
    this%astar=2.0d0*pi*this%astar
   
    !write lattice information
    write(6,'(/a)') 'real space lattice vectors (in nm)'
    write(6,'(3(3f12.6/))') transpose(this%areal)
    write(6,'(a)') 'reciprocal space lattice vectors (in nm-1)'
    write(6,'(3(3f12.6/))') transpose(this%astar)
    write(6,'(a)') 'Supercell real space lattice vectors (in nm)'
    write(6,'(3(3f12.6/))') transpose(this%a_sc)
    write(6,'(a/,2X,3L3/)') 'Supercell periodicity along each direction:',this%periodic

    !get sc_vec_period for quicker check of minimal distance between positions using supercell periodicity
    if(this%periodic(1))then
        allocate(i1,source=[-1,0,1])
    else
        allocate(i1,source=[0])
    endif
    if(this%periodic(2))then
        allocate(i2,source=[-1,0,1])
    else
        allocate(i2,source=[0])
    endif
    if(this%periodic(3))then
        allocate(i3,source=[-1,0,1])
    else
        allocate(i3,source=[0])
    endif
    allocate(this%sc_vec_period(3,size(i1)*size(i2)*size(i3)),source=0.0d0)   !real space vectors used to minimize distance using supercell periodicity
    i_vec=1
    do l=1,size(i3)
        tmp_vec(:,3)=this%a_sc(3,:)*i3(l)
        do j=1,size(i2)
            tmp_vec(:,2)=this%a_sc(2,:)*i2(j)
            do i=1,size(i1)
                tmp_vec(:,1)=this%a_sc(1,:)*i1(i)
                this%sc_vec_period(:,i_vec)=sum(tmp_vec,2)
                i_vec=i_vec+1
            enddo
        enddo
    enddo
end subroutine

subroutine init_order(this,cell,dim_modes)
    !initialize the lattice with the order parameter and cell properties
    class(lattice),intent(inout)    :: this
    type(t_cell), intent(in)        :: cell
    integer,intent(in)              :: dim_modes(number_different_order_parameters)

    this%cell=cell
    this%nmag=cell%num_mag()
    !obsolete orderparameter format
    this%dim_mode=sum(dim_modes)
    Call this%ordpar%init(this%dim_lat,this%nmag,this%dim_mode)
    
    !new orderparameter format
    if(dim_modes(1)>0) Call this%M%init(this%dim_lat,this%nmag,dim_modes(1))
    if(dim_modes(2)>0) Call this%E%init(this%dim_lat,this%nmag,dim_modes(2))
    if(dim_modes(3)>0) Call this%B%init(this%dim_lat,this%nmag,dim_modes(3))
    if(dim_modes(4)>0) Call this%T%init(this%dim_lat,this%nmag,dim_modes(4))
end subroutine

subroutine lattice_position(this,ind,pos)
    !slow routine to get position of individual index
    class(lattice)          :: this
    integer,intent(in)      :: ind(4)  !(i_x,i_y,i_y,i_at)
    real(8),intent(out)     :: pos(3)

    pos=matmul(real(ind(1:3)-1,8),this%areal)
    pos=pos+this%cell%atomic(ind(4))%position
end subroutine

subroutine min_diffvec_period(this,diff)
    !get minimal vector checking supercell perdiocity
    class(lattice),intent(in)   :: this
    real(8),intent(inout)       :: diff(3)

    real(8)     :: vec(3,size(this%sc_vec_period,2))
    real(8)     :: norm(size(this%sc_vec_period,2))
    integer     :: i

    do i=1,size(this%sc_vec_period,2)
        vec(:,i)=diff-this%sc_vec_period(:,i)
    enddo
    norm=norm2(vec,1)
    diff=vec(:,minloc(norm,1))
end subroutine

function lattice_position_difference(this,ind1,ind2)result(diff)
    !routine to get the difference vector between two individual atoms in index format
    !r_diff=pos_1-pos_2 (minimize with lattice vectors)
    class(lattice),intent(in)   :: this
    integer,intent(in)          :: ind1(4),ind2(4)  !(i_x,i_y,i_y,i_at)
    real(8)                     :: diff(3)
    real(8)                     :: pos1(3),pos2(3),norm

    Call this%pos_ind(ind1,pos1)
    Call this%pos_ind(ind2,pos2)
    diff=pos1-pos2
    norm=norm2(diff)
    if(norm>minval(norm2(this%a_sc,2)*0.5d0))then
    !check reduce length using symmetries 
        Call min_diffvec_period(this,diff)
    endif
end function


subroutine lattice_used_order(this,used)
    class(lattice),intent(in)   :: this
    logical,intent(out)         :: used(number_different_order_parameters)

    used=.False.
    if(this%M%dim_mode>0) used(1)=.True.
    if(this%E%dim_mode>0) used(2)=.True.
    if(this%B%dim_mode>0) used(3)=.True.
    if(this%T%dim_mode>0) used(4)=.True.
end subroutine

function index_m_1(this,indm)result(ind1)
    class(lattice),intent(in)   :: this
    integer,intent(in)          :: indm(:) 
    integer                     :: ind1
    integer     :: i

    ind1=indm(1)
    do i=2,size(indm)
        ind1=ind1+(indm(i)-1)*product(this%dim_lat(1:i-1))
    enddo
end function

function index_1_3(this,ind1)result(indm)
    class(lattice),intent(in)   :: this
    integer,intent(in)          :: ind1
    integer                     :: indm(3) 
    integer                     :: tmp,prod

    prod=product(this%dim_lat(:2))
    indm(3)=(ind1-1)/prod+1
    tmp=ind1-prod*(indm(3)-1)
    prod=this%dim_lat(1)
    indm(2)=(tmp-1)/prod+1
    indm(1)=(tmp-1)-prod*(indm(2)-1)+1
end function

function index_4_1(this,ind4)result(ind1)
    !get the correct in index from an (im,ix,iy,iz) array to the (3,Nmag*prod(Nlat)) magnetic array
    !order
    !MAKE SURE THAT THE 4-index array is in the correct order
    class(lattice),intent(in)   :: this
    integer,intent(in)          :: ind4(4) !(im,ix,iy,iz) 
    integer                     :: ind1
    integer     :: i
    integer     :: dim_full(4)

    dim_full(2:4)=this%dim_lat
    dim_full(1)=this%nmag
    ind1=ind4(1)
    do i=2,4
        ind1=ind1+(ind4(i)-1)*product(dim_full(1:i-1))
    enddo
end function

subroutine set_order_point(this,order,point)
    class(lattice),intent(in)   ::  this
    integer,intent(in)          ::  order
    real(8),pointer,intent(out) ::  point(:)

    select case(order)
    case(1)
        point=>this%M%all_modes
    case(2)
        point=>this%E%all_modes
    case(3)
        point=>this%B%all_modes
    case(4)
        point=>this%T%all_modes
    case default
        write(*,*) 'order:',order
        STOP 'failed to associate pointer in set_order_point'
    end select
end subroutine

subroutine set_order_point_single(this,order,ilat,point)
    class(lattice),intent(in)   :: this
    integer,intent(in)          :: order,ilat
    real(8),pointer,intent(out) :: point(:)

    integer                     :: bnd(2),dim_mode

    dim_mode=this%get_order_dim(order)
    bnd(1)=dim_mode*(ilat-1)+1
    bnd(2)=dim_mode*ilat
    select case(order)
    case(1)
        point=>this%M%all_modes(bnd(1):bnd(2))
    case(2)
        point=>this%E%all_modes(bnd(1):bnd(2))
    case(3)
        point=>this%B%all_modes(bnd(1):bnd(2))
    case(4)
        point=>this%T%all_modes(bnd(1):bnd(2))
    case default
        write(*,*) 'order:',order
        STOP 'failed to associate pointer in set_order_point'
    end select
end subroutine

subroutine reduce_cell(this,vec_in,order,vec_out)
    !adds all parts of a vectors parts of a vector together, that correspond to the same lattice cell
    class(lattice),intent(in)       ::  this
    real(8),intent(in)              ::  vec_in(:)
    integer,intent(in)              ::  order(:)
    real(8),intent(inout)           ::  vec_out(:)

    integer                 :: dim_modes(size(order))
    integer                 :: dim_mode_sum
    integer                 :: i

    do i=1,size(order)
        dim_modes(i)=this%get_order_dim(order(i))
    enddo
    dim_mode_sum=product(dim_modes)

    if(size(vec_in)/=this%Ncell*dim_mode_sum) STOP "vec_in of reduce has wrong dimension"
    if(size(vec_out)/=this%Ncell) STOP "vec_out of reduce has wrong dimension"

    vec_out=0.0d0
    do i=1,this%Ncell
       !probably not very efficient...
       vec_out(i)=sum(vec_in((i-1)*dim_mode_sum+1:i*dim_mode_sum))
    enddo
end subroutine

subroutine reduce(this,vec_in,order,order_keep,vec_out)
    class(lattice),intent(in)       ::  this
    real(8),intent(in)              ::  vec_in(:)
    integer,intent(in)              ::  order(:)
    integer,intent(in)              ::  order_keep
    real(8),intent(inout)           ::  vec_out(:)

    integer                 :: dim_modes(size(order))
    integer                 :: dim_mode_sum,dim_mode_keep
    integer                 :: ind_keep !index in dim_modes to keep (corresponding to index of points in set_order_comb)

    integer                 :: ind_div

    integer                 :: i
    integer                 :: site !actually site-1 for convenience for adding with site*dim_mode_sum
    integer                 :: dir

!   ind_keep=findloc(order,order_keep,dim=1)
    do i=1,size(order)
        if(order(i)==order_keep)then
            ind_keep=i
            exit
        endif
    enddo
    do i=1,size(order)
        dim_modes(i)=this%get_order_dim(order(i))
    enddo
    dim_mode_sum=product(dim_modes)
    dim_mode_keep=dim_modes(ind_keep)
    ind_div=product(dim_modes(:ind_keep-1))

    if(count(order==order_keep)/=1) STOP "implement reduce also for multiple entries of order_keep in order"
    if(size(vec_in)/=this%Ncell*dim_mode_sum) STOP "vec_in of reduce has wrong dimension"
    if(size(vec_out)/=this%Ncell*dim_modes(ind_keep)) STOP "vec_out of reduce has wrong dimension"

    vec_out=0.0d0
    do i=1,size(vec_in)
        site=(i-1)/dim_mode_sum
        dir=i-site*dim_mode_sum
        dir=modulo((dir-1)/ind_div,dim_mode_keep)+1
        vec_out(dir+site*dim_mode_keep)=vec_out(dir+site*dim_mode_keep)+vec_in(i)
    enddo
end subroutine

subroutine reduce_single(this,i_site,vec_in,order,order_keep,vec_out)
    class(lattice),intent(in)       :: this
    integer,intent(in)              :: i_site
    real(8),intent(in)              :: vec_in(:)
    integer,intent(in)              :: order(:)
    integer,intent(in)              :: order_keep
    real(8),intent(inout)           :: vec_out(:)

    integer                 :: dim_modes(size(order))
    integer                 :: dim_mode_sum,dim_mode_keep
    integer                 :: ind_keep !index in dim_modes to keep (corresponding to index of points in set_order_comb)

    integer                 :: ind_div

    integer                 :: i
    integer                 :: site !actually site-1 for convenience for adding with site*dim_mode_sum
    integer                 :: dir

!   ind_keep=findloc(order,order_keep,dim=1)
    do i=1,size(order)
        if(order(i)==order_keep)then
            ind_keep=i
            exit
        endif
    enddo
    do i=1,size(order)
        dim_modes(i)=this%get_order_dim(order(i))
    enddo
    dim_mode_sum=product(dim_modes)
    dim_mode_keep=dim_modes(ind_keep)
    ind_div=product(dim_modes(:ind_keep-1))

    if(count(order==order_keep)/=1) STOP "implement reduce also for multiple entries of order_keep in order"
    if(size(vec_in)/=dim_mode_sum) STOP "vec_in of reduce has wrong dimension"
    if(size(vec_out)/=dim_modes(ind_keep)) STOP "vec_out of reduce has wrong dimension"

    vec_out=0.0d0
    do i=1,size(vec_in)
        dir=modulo((i-1)/ind_div,dim_mode_keep)+1
        vec_out(dir)=vec_out(dir)+vec_in(i)
    enddo
end subroutine

subroutine set_order_comb_single(this,order,i_site,vec)
    !fills the combination of several order parameters according to order
    !probably quite slow implementation, but should at least work for any reasonable size of order
    !I SHOULD CHECK HOW SLOW THIS IS
    !vec should already be allocated to the size of the final vector ->product(dim_modes)
    use m_user_info
    class(lattice),intent(in)         ::  this
    integer,intent(in)                ::  order(:),i_site
    real(8),intent(inout)             ::  vec(:)

    type(point_arr)         :: points(size(order))
    integer                 :: dim_modes(size(order))
    integer                 :: dim_mode_sum
    integer                 :: i,i_entry,i_ord
    integer                 :: ind_site(size(order))
    integer                 :: ind(size(order))
    integer                 :: ind_div(size(order))
    integer                 :: entry_per_site
    
    do i=1,size(order)
        Call this%set_order_point(order(i),points(i)%v)
        dim_modes(i)=this%get_order_dim(order(i))
    enddo
    dim_mode_sum=product(dim_modes)
    do i=1,size(order)
        ind_div(i)=product(dim_modes(:i-1))
    enddo
    vec=1.0d0
    entry_per_site=product(dim_modes)
    ind_site=(i_site-1)*dim_modes
    do i=1,entry_per_site
        ind=(i-1)/ind_div
        ind=modulo(ind,dim_modes)+1+ind_site
        i_entry=i
        do i_ord=1,size(order)
            vec(i_entry)=vec(i_entry)*points(i_ord)%v(ind(i_ord))
        enddo
    enddo
end subroutine

subroutine set_order_comb(this,order,vec)
    !fills the combination of several order parameters according to order
    !probably quite slow implementation, but should at least work for any reasonable size of order
    !I SHOULD CHECK HOW SLOW THIS IS
    !vec should already be allocated to the size of the final vector ->dimH
    class(lattice),intent(in)         ::  this
    integer,intent(in)                ::  order(:)
    real(8),intent(inout)             ::  vec(:)
    !internal
    type(point_arr)         :: points(size(order))
    integer                 :: dim_modes(size(order))
    integer                 :: dim_mode_sum
    integer                 :: i_site,i,i_entry,i_ord
    integer                 :: ind_site(size(order))
    integer                 :: ind(size(order))
    integer                 :: ind_div(size(order))
    integer                 :: entry_per_site
    
    do i=1,size(order)
        Call this%set_order_point(order(i),points(i)%v)
        dim_modes(i)=this%get_order_dim(order(i))
    enddo
    dim_mode_sum=product(dim_modes)
    do i=1,size(order)
        ind_div(i)=product(dim_modes(:i-1))
    enddo
    vec=1.0d0
    entry_per_site=product(dim_modes)
    do i_site=1,this%ncell
        ind_site=(i_site-1)*dim_modes
        do i=1,entry_per_site
            ind=(i-1)/ind_div
            ind=modulo(ind,dim_modes)+1+ind_site
            i_entry=i+(i_site-1)*dim_mode_sum
            do i_ord=1,size(order)
                vec(i_entry)=vec(i_entry)*points(i_ord)%v(ind(i_ord))
            enddo
        enddo
    enddo
end subroutine

subroutine set_order_comb_exc(this,order,vec,order_exc)
    !fills the combination of several order parameters according to order
    !probably quite slow implementation, but should at least work for any reasonable size of order
    !TODO: I SHOULD CHECK HOW SLOW THIS IS
    !vec should already be allocated to the size of the final vector ->dimH
    class(lattice),intent(in)         ::  this
    integer,intent(in)                ::  order(:)
    logical,intent(in)                ::  order_exc(:)
    real(8),intent(inout)             ::  vec(:)

    type(point_arr)         :: points(size(order))
    integer                 :: dim_modes(size(order))
    integer                 :: dim_mode_sum
    integer                 :: i_site,i,i_entry,i_ord
    integer                 :: ind_site(size(order))
    integer                 :: ind(size(order))
    integer                 :: ind_div(size(order))

    if(size(order_exc)/=size(order)) STOP 'set_order_comb_exc input array shape wrong'
    do i=1,size(order)
        Call this%set_order_point(order(i),points(i)%v)
        dim_modes(i)=this%get_order_dim(order(i))
    enddo
    dim_mode_sum=product(dim_modes)
    do i=1,size(order)
        ind_div(i)=product(dim_modes(:i-1))
    enddo
    vec=1.0d0
    do i_site=1,this%ncell
        ind_site=(i_site-1)*dim_modes
        do i=1,product(dim_modes)
            ind=(i-1)/ind_div
            ind=modulo(ind,dim_modes)+1+ind_site
            i_entry=i+(i_site-1)*dim_mode_sum
            do i_ord=1,size(order)
                if(.not.order_exc(i_ord)) vec(i_entry)=vec(i_entry)*points(i_ord)%v(ind(i_ord)) 
            enddo
        enddo
    enddo
end subroutine

subroutine set_order_comb_exc_single(this,i_site,order,vec,order_exc)
    !fills the combination of several order parameters according to order
    !probably quite slow implementation, but should at least work for any reasonable size of order
    !TODO: I SHOULD CHECK HOW SLOW THIS IS
    !vec should already be allocated to the size of the final vector ->dimH
    class(lattice),intent(in)         ::  this
    integer,intent(in)                ::  order(:),i_site
    logical,intent(in)                ::  order_exc(:)
    real(8),intent(inout)             ::  vec(:)
    !internal
    type(point_arr)         :: points(size(order))
    integer                 :: dim_modes(size(order))
    integer                 :: dim_mode_sum
    integer                 :: i,i_entry,i_ord
    integer                 :: ind_site(size(order))
    integer                 :: ind(size(order))
    integer                 :: ind_div(size(order))

    if(size(order_exc)/=size(order)) STOP 'set_order_comb_exc input array shape wrong'
    do i=1,size(order)
        Call this%set_order_point(order(i),points(i)%v)
        dim_modes(i)=this%get_order_dim(order(i))
    enddo
    dim_mode_sum=product(dim_modes)
    do i=1,size(order)
        ind_div(i)=product(dim_modes(:i-1))
    enddo
    vec=1.0d0
    ind_site=(i_site-1)*dim_modes
    do i=1,product(dim_modes)
        ind=(i-1)/ind_div
        ind=modulo(ind,dim_modes)+1+ind_site
        i_entry=i
        do i_ord=1,size(order)
            if(.not.order_exc(i_ord)) vec(i_entry)=vec(i_entry)*points(i_ord)%v(ind(i_ord)) 
        enddo
    enddo
end subroutine

subroutine point_order_onsite(lat,op,dimH,modes,vec)
    !Subroutine that points modes to the order parameter vector according to op and dimH input
    !If size(op)>1 (i.e. dimension is folded from higher rank) allocates vec, sets it correctly
    !, and points modes
    !This only works if the unfolded order paramters are only considered on the same site
    class(lattice), intent(in)              :: lat
    integer,intent(in)                      :: op(:)
    integer,intent(in)                      :: dimH
    real(8),pointer,intent(out)             :: modes(:)
    real(8),allocatable,target,intent(out)  :: vec(:)

    if(size(op)==1)then
        Call lat%set_order_point(op(1),modes)
    else
        allocate(vec(dimH),source=0.0d0)
        Call lat%set_order_comb(op,vec)
        modes=>vec
    endif
end subroutine

subroutine point_order_onsite_single(lat,op,i_site,dim_mode,modes,vec)
    !Subroutine that points modes to the order parameter vector according to op and dimH input
    !If size(op)>1 (i.e. dimension is folded from higher rank) allocates vec, sets it correctly
    !, and points modes
    !This only works if the unfolded order paramters are only considered on the same site
    class(lattice), intent(in)              :: lat
    integer,intent(in)                      :: op(:)
    integer,intent(in)                      :: dim_mode,i_site
    real(8),pointer,intent(out)             :: modes(:)
    real(8),allocatable,target,intent(out)  :: vec(:)

    if(size(op)==1)then
        Call set_order_point_single(lat,op(1),i_site,modes)
    else
        allocate(vec(dim_mode),source=0.0d0)
        Call set_order_comb_single(lat,op,i_site,vec)
        modes=>vec
    endif
end subroutine

function get_order_dim(this,order) result(dim_mode)
    class(lattice),intent(in)   ::  this
    integer,intent(in)          ::  order
    integer                     ::  dim_mode

    select case(order)
    case(1)
        dim_mode=this%M%dim_mode
    case(2)
        dim_mode=this%E%dim_mode
    case(3)
        dim_mode=this%B%dim_mode
    case(4)
        dim_mode=this%T%dim_mode
    case default
        write(*,*) "order=",order
        STOP "trying to get dim_mode for unset orderparameter"
    end select
    if(dim_mode<1)then
        write(*,*) 'trying to get dim_mode for order=',order 
        ERROR STOP "dim_mode not positive, requested order parameter not set?"
    endif
end function

subroutine lattice_get_position_mag(this,pos)
    class(lattice),intent(in)       :: this
    real(8),allocatable,intent(out) :: pos(:,:,:,:,:)

    real(8), allocatable            :: pos_lat(:,:,:,:)
    integer, allocatable            :: ind_at_mag(:)
    integer     :: Nat_mag
    real(8)     :: pos_at(3) 
    integer     :: i,i1,i2,i3

    Call this%cell%ind_mag(ind_at_mag)
    Nat_mag=size(ind_at_mag)
    call get_position_cell(pos_lat,this%dim_lat,this%areal)
    allocate(pos(3,Nat_mag,this%dim_lat(1),this%dim_lat(2),this%dim_lat(3)),source=0.0d0)
    do i=1,Nat_mag  !could be done smarter with position array (:,Nat)
        pos_at=this%cell%atomic(ind_at_mag(i))%position
        pos(:,i,:,:,:)=pos_lat(:,:,:,:)
        do i3=1,this%dim_lat(3)
            do i2=1,this%dim_lat(2)
                do i1=1,this%dim_lat(1)
                    pos(:,i,i1,i2,i3)=pos(:,i,i1,i2,i3)+pos_at
                enddo
            enddo
        enddo
    enddo
    deallocate(pos_lat)
end subroutine

subroutine get_position_cell(pos,dim_lat,lat)
    !get the positions of the 0 point of each cell
    real(8), intent(out),allocatable    :: pos(:,:,:,:)
    integer, intent(in)                 :: dim_lat(:)
    real(8), intent(in)                 :: lat(:,:)
    ! internal variables
    real(8)     :: tmp(3,3)
    integer     :: i_z,i_y,i_x
    
    allocate(pos(3,dim_lat(1),dim_lat(2),dim_lat(3)),source=0.0d0)
    do i_z=1,dim_lat(3)
        tmp(:,3)=lat(3,:)*real(i_z-1,8)
        do i_y=1,dim_lat(2)
            tmp(:,2)=lat(2,:)*real(i_y-1,8)
            do i_x=1,dim_lat(1)
                tmp(:,1)=lat(1,:)*real(i_x-1,8)
                pos(1:3,i_x,i_y,i_z)=sum(tmp,dim=2)
            enddo
        enddo
    enddo
end subroutine 

subroutine copy_lattice(self,copy)
    class(lattice),intent(in)    :: self
    class(lattice),intent(out)   :: copy
    
    copy%areal=self%areal
    copy%astar=self%astar
    copy%a_sc=self%a_sc
    copy%dim_lat=self%dim_lat
    copy%ncell=self%ncell
    copy%n_system=self%n_system
    copy%dim_mode=self%dim_mode
    copy%nmag=self%nmag
    copy%periodic=self%periodic
    if(allocated(self%world)) allocate(copy%world,source=self%world)
    if(allocated(self%sc_vec_period)) allocate(copy%sc_vec_period,source=self%sc_vec_period)
    if(self%ordpar%dim_mode>0) Call self%ordpar%copy(copy%ordpar,self%dim_lat,self%nmag)

    if(self%m%dim_mode>0) Call self%M%copy(copy%M,self%dim_lat,self%nmag)
    if(self%e%dim_mode>0) Call self%E%copy(copy%E,self%dim_lat,self%nmag)
    if(self%b%dim_mode>0) Call self%B%copy(copy%B,self%dim_lat,self%nmag)
    if(self%t%dim_mode>0) Call self%T%copy(copy%T,self%dim_lat,self%nmag)
end subroutine

subroutine copy_val_lattice(self,copy)
    class(lattice),intent(inout)    :: self,copy
   
    !might make sense to check other values to make sure, but this should be faster
    Call self%ordpar%copy_val(copy%ordpar)

    if(self%m%dim_mode>0) Call self%M%copy_val(copy%M)
    if(self%e%dim_mode>0) Call self%E%copy_val(copy%E)
    if(self%b%dim_mode>0) Call self%B%copy_val(copy%B)
    if(self%t%dim_mode>0) Call self%T%copy_val(copy%T)
end subroutine

function op_name_to_int(name_in)result(int_out)
    character(len=*),intent(in) ::  name_in
    integer                     ::  int_out

    integer         ::  i
    
    do i=1,number_different_order_parameters
        if(trim(adjustl(name_in))==order_parameter_name(i))then
            int_out=i
            return
        endif
    enddo
    write(*,*) 'inserted name:',name_in
    STOP "Failed to identify name with operator"
end function

function op_int_to_abbrev(int_in)result(abbrev)
    integer,intent(in)          ::  int_in(:)   !operator integer description
    character(30)               ::  abbrev      

    integer         ::  i

    if(any(int_in>number_different_order_parameters)) ERROR STOP "too large int_in as operator description" 
    if(any(int_in<1)) ERROR STOP "too small int_in as operator description" 
    abbrev="" 
    do i=1,size(int_in)
        abbrev(i:i)=order_parameter_abbrev(int_in(i))
    enddo
end function
end module
