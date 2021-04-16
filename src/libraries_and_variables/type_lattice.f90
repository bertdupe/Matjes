module m_type_lattice
use m_basic_types
use m_order_parameter
use m_cell,only : t_cell
implicit none

private
public lattice, number_different_order_parameters,op_name_to_int,op_int_to_abbrev, op_abbrev_to_int,t_cell, order_par, order_parameter_name
public :: dim_modes_inner


integer,parameter :: number_different_order_parameters=5    !m,E,B,T,u
character(len=1),parameter :: order_parameter_abbrev(number_different_order_parameters)=["M","E","B","T","U"]
character(len=*),parameter :: order_parameter_name(number_different_order_parameters)=[&
                                &   'magnetic   ',&
                                &   'Efield     ',&
                                &   'Bfield     ',&
                                &   'temperature',&
                                &   'phonon     ']
integer,parameter           :: dim_modes_inner(number_different_order_parameters)=[3,3,3,1,3]   !rank of single mode entry (don't confuse with dim_modes which is super-cell resolved)

! variable that defines the lattice
type lattice
     real(8)        :: areal(3,3) !real space lattice vector
     real(8)        :: astar(3,3) !reciprocal space lattice vector
     type(t_cell)   :: cell
     integer        :: dim_lat(3) !number of unit-cells in each direction in supercell
     integer        :: Ncell !number of unit-cells (product of dim_lat)
     real(8)        :: a_sc(3,3) !real space lattice vectors of supercell
     real(8)        :: a_sc_inv(3,3) !inverse of real space lattice vectors of supercell (can be used to get back into supercell)
     integer        :: dim_modes(number_different_order_parameters)=0 !saves the dimension of the set order_parameters
     integer        :: site_per_cell(number_different_order_parameters)=0 !number of sites of order parameter per unit-cell (like nmag, nph)


     integer                :: n_system !no clue what this does
     integer, allocatable   :: world(:)
     logical                :: periodic(3) !lattice periodic in direction
     logical,private :: order_set(number_different_order_parameters)=.false. !logical variable which saves which orderparameters has been set

     !convenience parameters 
     integer        :: nmag=0 !number of magnetic atoms in the unit-cell (can also be obtained from cell)
     integer        :: nph=0  !number of atoms with phonons in the unit-cell (can also be obtained from cell) 
!     integer        :: nat=0  !overall number of atoms in the unit-cell (can also be obtained from cell) IMPLEMENT
!     integer        :: nattype=0  !overall number of different atom types (can also be obtained from cell) IMPLEMENT

     real(8),allocatable :: sc_vec_period(:,:)   !real space vectors used to minimize distance using supercell periodicity

!order parameters
     type(order_par)  :: M !magnetization 
     type(order_par)  :: E !Electric field
     type(order_par)  :: B !Magnetic field
     type(order_par)  :: T !Temperature
     type(order_par)  :: u !phonon

contains
    !initialization function
    procedure :: init_geo   !initialize the geometric properties (areal,Ncell,...)
    procedure :: init_order !initialize the order parameters dimensions
    procedure :: read_order !reads the order parameters from initialization files
    !basic function
    procedure :: copy => copy_lattice
    procedure :: copy_val_to => copy_val_lattice
    !mpi functions
    procedure :: bcast
    procedure :: bcast_val
    !get correct order parameter (or combination thereof)
    procedure :: set_order_point
    procedure :: set_order_comb
    procedure :: set_order_comb_exc
    procedure :: set_order_comb_exc_single
    procedure :: point_order => point_order_onsite
    procedure :: point_order_single => point_order_onsite_single
    procedure :: set_order_point_single_inner
    !!reduce that order parameter again
    procedure :: reduce
!    procedure :: reduce_single
    !!reduce vector in order parameter space to the site(e.g. energy Ncell resolution )
    procedure :: reduce_cell
    !real space position functions
    procedure :: get_pos_mag => lattice_get_position_mag
    procedure :: get_pos_center
    procedure :: pos_diff_ind => lattice_position_difference
    procedure :: pos_ind => lattice_position
    procedure :: dist_2vec  !distance between 2 vectors
    procedure :: min_diffvec  !get minimal difference vector using periodicity
    !index routines
    procedure :: index_m_1 
    procedure :: index_4_1 
    procedure :: index_1_3 
    procedure :: index_3_1_arr
    procedure :: index_4_1_arr 
    procedure :: fold_3_arr ! fold index 3 back to supercell
    !misc
    procedure :: used_order => lattice_used_order
    procedure :: get_order_dim 
end type lattice

type point_arr
    !private type used only here
    real(8),pointer     ::  v(:)
end type


contains 

subroutine init_geo(this,areal_in,alat,dim_lat,boundary)
    !initialize the lattice with the geometric inputs
    use m_vector, only : norm,cross
    use m_constants, only : pi
    class(lattice),intent(inout)    ::  this
    real(8),intent(in)              ::  areal_in(3,3),alat(3)
    integer,intent(in)              ::  dim_lat(3)
    logical,intent(in)              ::  boundary(3)
    !internal
    real(8)                         :: volume
    integer                         :: i
    !use supercell parameter
    integer,allocatable,dimension(:)     :: i1,i2,i3 
    real(8)         :: tmp_vec(3,3)
    integer         :: j,l,i_vec

    real(8)         ::  areal(3,3)

    do i=1,3
       this%areal(i,:)=areal_in(i,:)*alat(i)
    enddo
    areal=transpose(this%areal)
    this%dim_lat=dim_lat
    this%ncell=product(dim_lat)
    do i=1,3
       this%a_sc(i,:)=this%areal(i,:)*dim_lat(i)
    enddo

    this%periodic=boundary

    ! build up the reciprocal lattice vectors
    volume=dot_product(areal(:,1),cross(areal(:,2),areal(:,3)))
    this%astar(1,:) = cross(areal(:,2),areal(:,3))/volume
    this%astar(2,:) = cross(areal(:,3),areal(:,1))/volume
    this%astar(3,:) = cross(areal(:,1),areal(:,2))/volume
    this%astar=2.0d0*pi*this%astar

    this%a_sc_inv(1,:)= cross(this%a_sc(2,:),this%a_sc(3,:))
    this%a_sc_inv(2,:)= cross(this%a_sc(3,:),this%a_sc(1,:))
    this%a_sc_inv(3,:)= cross(this%a_sc(1,:),this%a_sc(2,:))
    this%a_sc_inv=this%a_sc_inv/dot_product(this%a_sc(1,:),cross(this%a_sc(2,:),this%a_sc(3,:)))   

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

subroutine init_order(this,cell,extpar_io)
    use m_input_types, only: extpar_input
    !initialize the lattice with the order parameter and cell properties
    class(lattice),intent(inout)    :: this
    type(t_cell), intent(in)        :: cell
    type(extpar_input),intent(in)   :: extpar_io

    integer                         :: i

    this%cell=cell
    this%nmag=count(cell%atomic(:)%moment/=0.0d0)
    this%nph=count(cell%atomic(:)%use_ph)
    !obsolete orderparameter format
    if(extpar_io%enable_M.and.this%nmag<1)then
        WRITE(*,'(///A)') 'CONGRATULATIONS, you have specified "enable_M" but set no magnetic atom'
        WRITE(*,'(A)') 'I though this should just do nothing and calculate as if there are no magnetic atoms'
        WRITE(*,'(A)') ', since there is nothing else this flag does other than throwing this error.'
        WRITE(*,'(A)') 'However Bertrand insists this should crash, so BAM'
        STOP 'Please tell Bertrand if you do not think this should crash'
    endif
    this%dim_modes(1)=this%nmag*3
    if(norm2(extpar_io%E)>0.0d0.or.extpar_io%enable_E) this%dim_modes(2)=3
    if(norm2(extpar_io%H)>0.0d0.or.extpar_io%enable_H) this%dim_modes(3)=3
    if(norm2(extpar_io%T)>0.0d0.or.extpar_io%enable_T) this%dim_modes(4)=1
    this%dim_modes(5)=this%nph*3
    !if(extpar_io%enable_u) dim_modes(5)=size(cell%atomic)*3
    this%order_set=this%dim_modes>0
    if(this%order_set(1)) this%site_per_cell(1)=this%nmag
    if(this%order_set(2)) this%site_per_cell(2)=1
    if(this%order_set(3)) this%site_per_cell(3)=1
    if(this%order_set(4)) this%site_per_cell(4)=1
    if(this%order_set(5)) this%site_per_cell(5)=this%nph

    !new orderparameter format
    if(this%order_set(1)) Call this%M%init(this%dim_lat,this%dim_modes(1),dim_modes_inner(1))
    if(this%order_set(2)) Call this%E%init(this%dim_lat,this%dim_modes(2),dim_modes_inner(2))
    if(this%order_set(3)) Call this%B%init(this%dim_lat,this%dim_modes(3),dim_modes_inner(3))
    if(this%order_set(4)) Call this%T%init(this%dim_lat,this%dim_modes(4),dim_modes_inner(4))
    if(this%order_set(5)) Call this%u%init(this%dim_lat,this%dim_modes(5),dim_modes_inner(5))
end subroutine

subroutine read_order(this,suffix_in,isinit_opt)
    !read the order parameters from a 
    class(lattice),intent(inout)        :: this
    character(*),intent(in),optional    :: suffix_in    !suffix to order parameter name for input file
    logical,intent(inout),optional      :: isinit_opt(number_different_order_parameters)    !(in) true is already initialized, (out) true if previously or now initialized

    logical                             :: isinit(number_different_order_parameters)
    character(*),parameter              :: suffix_default='.inp'
    character(:), allocatable           :: suffix

    if(present(suffix_in))then
        suffix=trim(adjustl(suffix_in))
    else
        suffix=suffix_default
    endif
    isinit=.false.
    if(present(isinit_opt)) isinit=isinit_opt
    if(this%order_set(1).and..not.isinit(1)) Call this%M%read_file(trim(order_parameter_name(1))//suffix,isinit(1))
    if(this%order_set(2).and..not.isinit(2)) Call this%E%read_file(trim(order_parameter_name(2))//suffix,isinit(2))
    if(this%order_set(3).and..not.isinit(3)) Call this%B%read_file(trim(order_parameter_name(3))//suffix,isinit(3))
    if(this%order_set(4).and..not.isinit(4)) Call this%T%read_file(trim(order_parameter_name(4))//suffix,isinit(4))
    if(this%order_set(5).and..not.isinit(5)) Call this%u%read_file(trim(order_parameter_name(5))//suffix,isinit(5))
    deallocate(suffix)
    if(present(isinit_opt)) isinit_opt=isinit
end subroutine

subroutine lattice_position(this,ind,pos)
    !slow routine to get position of individual index
    class(lattice)          :: this
    integer,intent(in)      :: ind(4)  !(i_x,i_y,i_y,i_at)
    real(8),intent(out)     :: pos(3)

    pos=matmul(real(ind(1:3)-1,8),this%areal)
    pos=pos+this%cell%atomic(ind(4))%position
end subroutine


subroutine dist_2vec(this,vec1,vec2,dist)
    !get the minimial distance between 2 vectors on the lattice obeying the symmetries
    class(lattice),intent(in)   :: this
    real(8),intent(in)          :: vec1(3),vec2(3)
    real(8),intent(out)         :: dist 
    !internal
    real(8)     :: vec(3,size(this%sc_vec_period,2))
    real(8)     :: norm(size(this%sc_vec_period,2))
    real(8)     :: vec1_sc(3),vec2_sc(3)
    integer     :: i

    vec1_sc=pos_first_sc(this,vec1)
    vec2_sc=pos_first_sc(this,vec2)
    do i=1,size(this%sc_vec_period,2)
        vec(:,i)=vec2_sc-vec1_sc-this%sc_vec_period(:,i)
    enddo
    norm=norm2(vec,1)
    dist=minval(norm)
end subroutine

subroutine min_diffvec(this,diff)
    !get minimal vector checking supercell perdiocity
    class(lattice),intent(in)   :: this
    real(8),intent(inout)       :: diff(3)

    real(8)     :: diff_sc(3)
    real(8)     :: vec(3,size(this%sc_vec_period,2))
    real(8)     :: norm(size(this%sc_vec_period,2))
    integer     :: i

    diff_sc=pos_first_sc(this,diff)
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
        Call this%min_diffvec(diff)
    endif
end function


subroutine lattice_used_order(this,used)
    class(lattice),intent(in)   :: this
    logical,intent(out)         :: used(number_different_order_parameters)

    used=this%order_set
end subroutine

subroutine index_3_1_arr(this,N,ind3,ind1)
    class(lattice),intent(in)   :: this
    integer,intent(in)          :: N
    integer,intent(inout)       :: ind3(3,N) !doesn't actually change, but gets used in between
    integer,intent(inout)       :: ind1(N)
    integer     :: i
    integer     :: mat(3)

    ind3=ind3-1
    mat=[1,this%dim_lat(1),this%dim_lat(1)*this%dim_lat(2)]
    ind1=matmul(mat,ind3)+1
    ind3=ind3+1
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

subroutine index_4_1_arr(this,N,ind4,ind1)
    !get the correct in index from an (im,ix,iy,iz) array to the (3,Nmag*prod(Nlat)) magnetic array
    !order
    !MAKE SURE THAT THE 4-index array is in the correct order
    class(lattice),intent(in)   :: this
    integer,intent(in)          :: N
    integer,intent(in)          :: ind4(4,N) !(im,ix,iy,iz) 
    integer,intent(out)         :: ind1(N)
    integer     :: i
    integer     :: dim_full(4)

    dim_full(2:4)=this%dim_lat
    dim_full(1)=this%nmag
    ind1=ind4(1,:)
    do i=2,4
        ind1=ind1+(ind4(i,:)-1)*product(dim_full(1:i-1))
    enddo
end subroutine 

subroutine fold_3_arr(this,N,ind_arr)
    !fold the lattice indices back, assuming the periodicity exists...
    class(lattice),intent(in)   ::  this
    integer,intent(in)          ::  N
    integer,intent(inout)       ::  ind_arr(3,N)
    integer                     ::  dim_lat_spread(3,N)

    dim_lat_spread=spread(this%dim_lat,2,N)

    ind_arr=ind_arr-((ind_arr-dim_lat_spread)/dim_lat_spread)*dim_lat_spread   !return periodically to positive entries
    ind_arr=ind_arr-((ind_arr-1)/dim_lat_spread)*dim_lat_spread                !return periodically to (1:dim_lat)
end subroutine


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
    case(5)
        point=>this%u%all_modes
    case default
        write(*,*) 'order:',order
        STOP 'failed to associate pointer in set_order_point'
    end select
end subroutine

subroutine set_order_point_single_inner(this,order,i_inner,point,bnd)
    !sets the pointer point to the 
    class(lattice),intent(in)   :: this
    integer,intent(in)          :: order,i_inner
    real(8),pointer,intent(out) :: point(:)
    integer,intent(out)         :: bnd(2)

    bnd=[1+(i_inner-1)*dim_modes_inner(order),i_inner*dim_modes_inner(order)]

    select case(order)
    case(1)
        point=>this%M%all_modes(bnd(1):bnd(2))
    case(2)
        point=>this%E%all_modes(bnd(1):bnd(2))
    case(3)
        point=>this%B%all_modes(bnd(1):bnd(2))
    case(4)
        point=>this%T%all_modes(bnd(1):bnd(2))
    case(5)
        point=>this%u%all_modes(bnd(1):bnd(2))
    case default
        write(*,*) 'order:',order
        STOP 'failed to associate pointer in set_order_point'
    end select
end subroutine

subroutine set_order_point_single(this,order,i_site,dim_bnd,dim_mode,point,bnd)
    class(lattice),intent(in)   :: this
    integer,intent(in)          :: order,i_site,dim_mode
    integer, intent(in)         :: dim_bnd(2,number_different_order_parameters)  
    real(8),pointer,intent(out) :: point(:)
    integer,intent(out)         :: bnd(2)

    Call get_bnd_single_1(this,order,i_site,dim_bnd,bnd)
    !this select case might be ultra slow, make another indexable pointer array in lattice?
    select case(order)
    case(1)
        point=>this%M%all_modes(bnd(1):bnd(2))
    case(2)
        point=>this%E%all_modes(bnd(1):bnd(2))
    case(3)
        point=>this%B%all_modes(bnd(1):bnd(2))
    case(4)
        point=>this%T%all_modes(bnd(1):bnd(2))
    case(5)
        point=>this%u%all_modes(bnd(1):bnd(2))
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

    integer                 :: dim_mode_sum
    integer                 :: i

    dim_mode_sum=product(this%dim_modes(order))

    if(size(vec_in)/=this%Ncell*dim_mode_sum) ERROR STOP "vec_in of reduce has wrong dimension"
    if(size(vec_out)/=this%Ncell) ERROR STOP "vec_out of reduce has wrong dimension"

    vec_out=0.0d0
    do i=1,this%Ncell
       !probably not very efficient...
       vec_out(i)=sum(vec_in((i-1)*dim_mode_sum+1:i*dim_mode_sum))
    enddo
end subroutine

subroutine reduce(this,vec_in,order,ind_keep,vec_out)
    class(lattice),intent(in)       :: this
    real(8),intent(in)              :: vec_in(:)
    integer,intent(in)              :: order(:)
    integer,intent(in)              :: ind_keep !index in dim_modes to keep (corresponding to index of points in set_order_comb)
    real(8),intent(inout)           :: vec_out(:)

    integer                 :: dim_mode_sum,dim_mode_keep
    integer                 :: dim_modes(size(order))   !dimension of modes supplied in order

    integer                 :: ind_div
    integer                 :: i
    integer                 :: site !actually site-1 for convenience for adding with site*dim_mode_sum
    integer                 :: dir

    dim_modes=this%dim_modes(order)
    dim_mode_sum=product(dim_modes)
    dim_mode_keep=dim_modes(ind_keep)
    ind_div=product(dim_modes(:ind_keep-1))

    if(size(vec_in)/=this%Ncell*dim_mode_sum) ERROR STOP "vec_in of reduce has wrong dimension"
    if(size(vec_out)/=this%Ncell*dim_modes(ind_keep)) STOP "vec_out of reduce has wrong dimension"

    vec_out=0.0d0
    do i=1,size(vec_in)
        site=(i-1)/dim_mode_sum
        dir=i-site*dim_mode_sum
        dir=modulo((dir-1)/ind_div,dim_mode_keep)+1
        vec_out(dir+site*dim_mode_keep)=vec_out(dir+site*dim_mode_keep)+vec_in(i)
    enddo
end subroutine

!subroutine reduce_single(this,i_site,vec_in,order,order_keep,vec_out)
!    class(lattice),intent(in)       :: this
!    integer,intent(in)              :: i_site
!    real(8),intent(in)              :: vec_in(:)
!    integer,intent(in)              :: order(:)
!    integer,intent(in)              :: order_keep
!    real(8),intent(inout)           :: vec_out(:)
!
!    integer                 :: dim_mode_sum,dim_mode_keep
!    integer                 :: ind_keep !index in dim_modes to keep (corresponding to index of points in set_order_comb)
!    integer                 :: dim_modes(size(order))   !dimension of modes supplied in order
!
!    integer                 :: ind_div
!
!    integer                 :: i
!    integer                 :: site !actually site-1 for convenience for adding with site*dim_mode_sum
!    integer                 :: dir
!
!!   ind_keep=findloc(order,order_keep,dim=1)
!    do i=1,size(order)
!        if(order(i)==order_keep)then
!            ind_keep=i
!            exit
!        endif
!    enddo
!    dim_modes=this%dim_modes(order)
!    dim_mode_sum=product(dim_modes)
!    dim_mode_keep=dim_modes(ind_keep)
!    ind_div=product(dim_modes(:ind_keep-1))
!
!    if(count(order==order_keep)/=1) ERROR STOP "implement reduce also for multiple entries of order_keep in order"
!    if(size(vec_in)/=dim_mode_sum) Error STOP "vec_in of reduce has wrong dimension"
!    if(size(vec_out)/=dim_modes(ind_keep)) ERROR STOP "vec_out of reduce has wrong dimension"
!
!    vec_out=0.0d0
!    do i=1,size(vec_in)
!        dir=modulo((i-1)/ind_div,dim_mode_keep)+1
!        vec_out(dir)=vec_out(dir)+vec_in(i)
!    enddo
!end subroutine

subroutine set_order_comb_single(this,order,i_site,dim_bnd_in,vec,bnd)
    !fills the combination of several order parameters according to order
    !probably quite slow implementation, but should at least work for any reasonable size of order
    !I SHOULD CHECK HOW SLOW THIS IS
    !vec should already be allocated to the size of the final vector ->product(dim_modes)
    class(lattice),intent(in)   :: this
    integer,intent(in)          :: order(:) !order integers (1:rank reduce operator)->(1:number_different_order_parameters)
    integer,intent(in)          :: i_site   !site index (1:Ncell)
    integer, intent(in)         :: dim_bnd_in(2,number_different_order_parameters)  
    real(8),intent(inout)       :: vec(:)
    integer,intent(out)         :: bnd(2)

    !local instance of dim_bnd to only have reduce bnd for the first occurance of the order-parameter (multiple M on one site?, biquadratic?)
    !(NOT TESTED SINCE NO SUCH HAMILTIONIAN EXISTS YES)
    integer                 :: dim_bnd(2,number_different_order_parameters)
    type(point_arr)         :: points(size(order))      !pointers to the respective order parameter
    integer                 :: dim_modes(size(order))   !overall dim_mode  (1:size(order))-> dim_mode of this order parameter
    integer                 :: dim_mode_all             !product of all considered order parameters
    integer                 :: ind_site(size(order))    !index offset to entries of sites for all considered order-parameters
    integer                 :: ind(size(order))         !index of considered component in space of respective order parameter
    integer                 :: ind_div(size(order))     !temp. values to get all indices as vector operation 
    integer                 :: bnd_loc(2,size(order))   !bondaries locally used for each 

    integer                 :: i,i_ord
    
    dim_bnd=dim_bnd_in
    do i=1,size(order)
        Call this%set_order_point(order(i),points(i)%v)
        dim_modes(i)=this%get_order_dim(order(i))
        bnd_loc(:,i)=dim_bnd(:,order(i))
        dim_bnd(:,order(i))=[1,dim_modes(i)]
    enddo
    dim_mode_all=product(dim_modes)
    do i=1,size(order)
        ind_div(i)=product(dim_modes(:i-1))
    enddo
    vec=1.0d0
    ind_site=(i_site-1)*dim_modes
    do i=1,dim_mode_all
        ind=(i-1)/ind_div
        ind=modulo(ind,dim_modes)+1
        if(any(bnd_loc(1,:)>ind).or.any(bnd_loc(2,:)<ind))then   !probably rather slow and certainly could be solved quicker
            vec(i)=0.0d0
            cycle
        endif
        ind=ind+ind_site
        do i_ord=1,size(order)
            vec(i)=vec(i)*points(i_ord)%v(ind(i_ord))
        enddo
    enddo
    bnd=(i_site-1)*dim_mode_all+[1,dim_mode_all]    !return full array, could be optimizes especially if the boundary order parameter is the slowly varying parameter
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
    integer                 :: dim_mode_all
    integer                 :: i_site,i,i_entry,i_ord
    integer                 :: ind_site(size(order))
    integer                 :: ind(size(order))
    integer                 :: ind_div(size(order))
    
    do i=1,size(order)
        Call this%set_order_point(order(i),points(i)%v)
        dim_modes(i)=this%get_order_dim(order(i))
    enddo
    dim_mode_all=product(dim_modes)
    do i=1,size(order)
        ind_div(i)=product(dim_modes(:i-1))
    enddo
    vec=1.0d0
    do i_site=1,this%ncell
        ind_site=(i_site-1)*dim_modes
        do i=1,dim_mode_all
            ind=(i-1)/ind_div
            ind=modulo(ind,dim_modes)+1+ind_site
            i_entry=i+(i_site-1)*dim_mode_all
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

subroutine point_order_onsite_single(lat,op,i_site,dim_bnd,dim_mode,modes,vec,bnd)
    !Subroutine that points modes to the order parameter vector according to op and dimH input
    !If size(op)>1 (i.e. dimension is folded from higher rank) allocates vec, sets it correctly
    !, and points modes
    !This only works if the unfolded order paramters are only considered on the same site
    class(lattice), intent(in)              :: lat
    integer,intent(in)                      :: op(:)
    integer,intent(in)                      :: dim_mode,i_site
    integer, intent(in)                     :: dim_bnd(2,number_different_order_parameters)    !starting/final index in respective dim_mode of the order parameter (so that energy of single magnetic atom can be be calculated
    real(8),pointer,intent(out)             :: modes(:)
    real(8),allocatable,target,intent(out)  :: vec(:)   !space to allocate array if not single operator
    integer,intent(out)                     :: bnd(2)   !boundary indices in (1:dim_mode*N_cell)-basis

    if(size(op)==1)then
        Call set_order_point_single(lat,op(1),i_site,dim_bnd,dim_mode,modes,bnd)
    else
        allocate(vec(dim_mode),source=0.0d0)
        Call set_order_comb_single(lat,op,i_site,dim_bnd,vec,bnd)
        modes=>vec
    endif
end subroutine

subroutine get_bnd_single_1(lat,op,i_site,dim_bnd,bnd)
    class(lattice), intent(in)  :: lat
    integer,intent(in)          :: op
    integer,intent(in)          :: i_site
    integer,intent(in)          :: dim_bnd(2,number_different_order_parameters)
    integer,intent(out)         :: bnd(2)

    bnd=dim_bnd(:,op)+(i_site-1)*lat%dim_modes(op)
end subroutine

subroutine get_bnd_single(lat,op,i_site,dim_bnd,bnd)
    class(lattice), intent(in)   :: lat
    integer,intent(in)          :: op(:)
    integer,intent(in)          :: i_site
    integer,intent(in)          :: dim_bnd(2,number_different_order_parameters)
    integer,intent(out)         :: bnd(2)

    integer     ::  dim_mode

    if(size(op)==1)then
        Call get_bnd_single_1(lat,op(1),i_site,dim_bnd,bnd)
    else
        STOP "IMPLEMENT FOR RANK2 and higher parameters"
    endif

end subroutine


function get_order_dim(this,order,ignore) result(dim_mode)
    class(lattice),intent(in)   :: this
    integer,intent(in)          :: order
    logical,intent(in),optional :: ignore
    integer                     :: dim_mode

#ifdef CPP_DEBUG
    !this could spot a programming mistake
    if(order<1.or.order>number_different_order_parameters)then
        write(*,*) "order=",order
        STOP "trying to get dim_mode for unset orderparameter"
    endif
#endif
    dim_mode=this%dim_modes(order)
    if(dim_mode<1)then
        if(present(ignore))then
            if(ignore) return
        endif
        write(*,*) 'trying to get dim_mode for order=',order 
        ERROR STOP "dim_mode not positive, requested order parameter not set?"
    endif
end function

subroutine get_pos_center(this,pos_out)
    class(lattice),intent(in)               :: this
    real(8),allocatable,intent(out),target  :: pos_out(:)

    real(8),pointer,contiguous      :: pos(:,:,:,:)
    real(8),pointer,contiguous      :: pos3(:,:)
    real(8),allocatable             :: pos_lat(:,:,:,:)
    real(8)                         :: pos_center(3) 
    integer                         :: i

    call get_position_cell(pos_lat,this%dim_lat,this%areal)
    allocate(pos_out(3*product(this%dim_lat)),source=0.0d0)
    pos(1:3,1:this%dim_lat(1),1:this%dim_lat(2),1:this%dim_lat(3))=>pos_out
    pos3(1:3,1:product(this%dim_lat))=>pos_out
    pos_center=sum(this%areal,1)*0.5d0
    pos=pos_lat
    do i=1,size(pos3,2)
        pos3(:,i)=pos3(:,i)+pos_center
    enddo
    deallocate(pos_lat)
    nullify(pos,pos3)
end subroutine

subroutine lattice_get_position_mag(this,pos_out)
    class(lattice),intent(in)               :: this
    real(8),allocatable,intent(out),target  :: pos_out(:)

    real(8),pointer                 :: pos(:,:,:,:,:)
    real(8), allocatable            :: pos_lat(:,:,:,:)
    integer, allocatable            :: ind_at_mag(:)
    integer     :: Nat_mag
    real(8)     :: pos_at(3) 
    integer     :: i,i1,i2,i3

    Call this%cell%ind_mag_all(ind_at_mag)
    Nat_mag=size(ind_at_mag)
    call get_position_cell(pos_lat,this%dim_lat,this%areal)
    allocate(pos_out(3*Nat_mag*product(this%dim_lat)),source=0.0d0)
    pos(1:3,1:Nat_mag,1:this%dim_lat(1),1:this%dim_lat(2),1:this%dim_lat(3))=>pos_out
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
    nullify(pos)
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

subroutine copy_lattice(this,copy)
    class(lattice),intent(in)    :: this
    class(lattice),intent(out)   :: copy
    
    copy%areal=this%areal
    copy%astar=this%astar
    copy%a_sc=this%a_sc
    copy%a_sc_inv=this%a_sc_inv
    copy%dim_lat=this%dim_lat
    copy%ncell=this%ncell
    copy%n_system=this%n_system
    copy%dim_modes=this%dim_modes
    copy%nmag=this%nmag
    copy%nph=this%nph
    copy%periodic=this%periodic
    if(allocated(this%world)) allocate(copy%world,source=this%world)
    if(allocated(this%sc_vec_period)) allocate(copy%sc_vec_period,source=this%sc_vec_period)
    copy%order_set=this%order_set
    copy%site_per_cell=this%site_per_cell
    Call this%cell%copy(copy%cell)

    if(this%order_set(1)) Call this%M%copy(copy%M,this%dim_lat)
    if(this%order_set(2)) Call this%E%copy(copy%E,this%dim_lat)
    if(this%order_set(3)) Call this%B%copy(copy%B,this%dim_lat)
    if(this%order_set(4)) Call this%T%copy(copy%T,this%dim_lat)
    if(this%order_set(5)) Call this%u%copy(copy%u,this%dim_lat)
end subroutine

subroutine copy_val_lattice(this,copy)
    class(lattice),intent(inout)    :: this,copy
   
    if(this%order_set(1)) Call this%M%copy_val(copy%M)
    if(this%order_set(2)) Call this%E%copy_val(copy%E)
    if(this%order_set(3)) Call this%B%copy_val(copy%B)
    if(this%order_set(4)) Call this%T%copy_val(copy%T)
    if(this%order_set(5)) Call this%u%copy_val(copy%u)
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

function op_abbrev_to_int(abbrev_in)result(int_out)
    character(len=*), intent(in)    :: abbrev_in
    integer                         :: int_out(len(abbrev_in))   !operator integer description
    integer         ::  i,j

    outer:do i=1,len(abbrev_in)
        do j=1,size(order_parameter_abbrev)
            if(abbrev_in(i:i)==order_parameter_abbrev(j))then
                int_out(i)=j
                cycle outer
            endif
        enddo
        write(*,*) "failing to get operator integer for ",abbrev_in(i:i)
        STOP "DID NOT FIND ABBREVIATION"
    enddo outer
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

function pos_first_sc(this,vec_in)result(vec_out)
    !returns the vector position back in the first supercell
    class(lattice),intent(in)   :: this
    real(8),intent(in)          :: vec_in(3)
    real(8)                     :: vec_out(3)
    integer                     :: ind(3)

    ind=floor(matmul(this%a_sc_inv,vec_in))
    vec_out=vec_in-matmul(real(ind,8),this%a_sc)
end function

subroutine bcast(this,comm)
use mpi_basic                
    class(lattice),intent(inout)    ::  this
    type(mpi_type),intent(in)       ::  comm

#ifdef CPP_MPI
    integer     :: ierr
    integer     :: N

    CALL this%cell%bcast(comm)
    Call MPI_Bcast(this%areal    , 9                   , MPI_REAL8  , comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%astar    , 9                   , MPI_REAL8  , comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%dim_lat  , 3                   , MPI_INTEGER, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%Ncell    , 1                   , MPI_INTEGER, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%a_sc     , 9                   , MPI_REAL8  , comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%a_sc_inv , 9                   , MPI_REAL8  , comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%dim_modes, size(this%dim_modes), MPI_INTEGER, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%n_system , 1                   , MPI_INTEGER, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%nmag     , 1                   , MPI_INTEGER, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%nph      , 1                   , MPI_INTEGER, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%periodic , 3                   , MPI_LOGICAL, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%order_set, size(this%order_set), MPI_LOGICAL, comm%mas, comm%com,ierr)

    if(comm%ismas) N=size(this%world)
    Call MPI_Bcast(N, 1, MPI_INTEGER, comm%mas, comm%com,ierr)
    if(.not.comm%ismas) allocate(this%world(N))
    Call MPI_Bcast(this%world, N, MPI_INTEGER, comm%mas, comm%com,ierr)

    if(comm%ismas) N=size(this%world)
    Call MPI_Bcast(N, 1, MPI_INTEGER, comm%mas, comm%com,ierr)
    if(.not.comm%ismas) allocate(this%sc_vec_period(3,N/3))
    Call MPI_Bcast(this%sc_vec_period, N, MPI_REAL8, comm%mas, comm%com,ierr)

    if(this%order_set(1)) Call this%M%bcast(comm)
    if(this%order_set(2)) Call this%E%bcast(comm)
    if(this%order_set(3)) Call this%B%bcast(comm)
    if(this%order_set(4)) Call this%T%bcast(comm)
    if(this%order_set(5)) Call this%u%bcast(comm)
#else
    continue
#endif
end subroutine

subroutine bcast_val(this,comm)
use mpi_basic                
    class(lattice),intent(inout)    ::  this
    type(mpi_type),intent(in)       ::  comm

#ifdef CPP_MPI
    integer     :: ierr
    integer     :: N

    if(this%order_set(1)) Call this%M%bcast_val(comm)
    if(this%order_set(2)) Call this%E%bcast_val(comm)
    if(this%order_set(3)) Call this%B%bcast_val(comm)
    if(this%order_set(4)) Call this%T%bcast_val(comm)
    if(this%order_set(5)) Call this%u%bcast_val(comm)
#else
    continue
#endif
end subroutine

end module
