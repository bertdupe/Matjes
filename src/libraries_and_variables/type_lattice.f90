module m_type_lattice
use m_basic_types
use m_order_parameter
use m_cell,only : t_cell
use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit, output_unit
implicit none

private
public lattice, number_different_order_parameters,op_name_to_int,op_int_to_abbrev, op_abbrev_to_int,t_cell, order_par, order_parameter_name
public :: dim_modes_inner


integer,parameter :: number_different_order_parameters=6    !m,E,B,T,u,w
character(len=1),parameter :: order_parameter_abbrev(number_different_order_parameters)=["M","E","B","T","U","W"]
character(len=*),parameter :: order_parameter_name(number_different_order_parameters)=[&
                                &   'magnetic    ',&
                                &   'Efield      ',&
                                &   'Bfield      ',&
                                &   'temperature ',&
                                &   'phonon      ',&
                                &   'wavefunction']

integer,parameter           :: dim_modes_inner(number_different_order_parameters)=[3,3,3,1,3,2]   !rank of single mode entry (don't confuse with dim_modes which is super-cell resolved)

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
     integer        :: norb=0  !number of wave function orbitals per cell 
!     integer        :: nat=0  !overall number of atoms in the unit-cell (can also be obtained from cell) IMPLEMENT
!     integer        :: nattype=0  !overall number of different atom types (can also be obtained from cell) IMPLEMENT

     real(8),allocatable :: sc_vec_period(:,:)   !real space vectors used to minimize distance using supercell periodicity

!order parameters
     type(order_par)  :: M !magnetization 
     type(order_par)  :: E !Electric field
     type(order_par)  :: B !Magnetic field
     type(order_par)  :: T !Temperature
     type(order_par)  :: u !phonon
     type(order_par)  :: w !wavefunction

contains
    !initialization function
    procedure :: init_geo   !initialize the geometric properties (areal,Ncell,...)
    procedure :: init_order !initialize the order parameters dimensions
    procedure :: read_order !reads the order parameters from initialization files
    !basic function
    procedure :: copy => copy_lattice
    procedure :: copy_val_to => copy_val_lattice
    !netcdf io functions
    procedure :: io_open
    procedure :: io_write
    procedure :: io_close
    !mpi functions
    procedure :: bcast
    procedure :: bcast_val
    !get correct order parameter (or combination thereof)
    procedure :: set_order_point
    procedure :: set_order_point_inner
!    procedure :: get_pos_mag => lattice_get_position_mag
!    procedure :: get_pos_center
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
    real(8),pointer,contiguous     ::  v(:)
end type


contains 

subroutine io_open(this)
    class(lattice),intent(inout)        ::  this
    integer ::  i,ind

    if(this%order_set(1)) call this%m%open_io(trim(order_parameter_name(1)))
    if(this%order_set(2)) call this%e%open_io(trim(order_parameter_name(2)))
    if(this%order_set(3)) call this%b%open_io(trim(order_parameter_name(3)))
    if(this%order_set(4)) call this%t%open_io(trim(order_parameter_name(4)))
    if(this%order_set(5)) call this%u%open_io(trim(order_parameter_name(5)))
    if(this%order_set(5)) call this%w%open_io(trim(order_parameter_name(5)))
end subroutine

subroutine io_write(this)
    class(lattice),intent(inout)        ::  this
    integer ::  i,ind

    if(this%order_set(1)) call this%m%write_io()
    if(this%order_set(2)) call this%e%write_io()
    if(this%order_set(3)) call this%b%write_io()
    if(this%order_set(4)) call this%t%write_io()
    if(this%order_set(5)) call this%u%write_io()
    if(this%order_set(6)) call this%w%write_io()
end subroutine

subroutine io_close(this)
    class(lattice),intent(inout)        ::  this
    integer ::  i,ind

    if(this%order_set(1)) call this%m%close_io()
    if(this%order_set(2)) call this%e%close_io()
    if(this%order_set(3)) call this%b%close_io()
    if(this%order_set(4)) call this%t%close_io()
    if(this%order_set(5)) call this%u%close_io()
    if(this%order_set(6)) call this%w%close_io()
end subroutine


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
    real(8)         ::  a_sc(3,3)

    areal=transpose(areal_in)
    do i=1,3
       areal(:,i)=areal(:,i)*alat(i)
    enddo
    this%areal=transpose(areal)

    this%dim_lat=dim_lat
    this%ncell=product(dim_lat)
    do i=1,3
       a_sc(:,i)=areal(:,i)*dim_lat(i)
    enddo
    this%a_sc=transpose(a_sc)

    this%periodic=boundary

    ! build up the reciprocal lattice vectors
    volume=dot_product(areal(:,1),cross(areal(:,2),areal(:,3)))
    this%astar(1,:) = cross(areal(:,2),areal(:,3))/volume
    this%astar(2,:) = cross(areal(:,3),areal(:,1))/volume
    this%astar(3,:) = cross(areal(:,1),areal(:,2))/volume
    this%astar=2.0d0*pi*this%astar

    this%a_sc_inv(1,:)= cross(a_sc(:,2),a_sc(:,3))
    this%a_sc_inv(2,:)= cross(a_sc(:,3),a_sc(:,1))
    this%a_sc_inv(3,:)= cross(a_sc(:,1),a_sc(:,2))
    this%a_sc_inv=this%a_sc_inv/dot_product(a_sc(:,1),cross(a_sc(:,2),a_sc(:,3)))   

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
    this%norb=sum(cell%atomic(:)%orbitals)
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
    if(extpar_io%enable_w)then
        this%dim_modes(6)=this%norb*2    !requires to be enables manually because of separate TB-part no requiring the lattice, factor 2 because is complex
        if(any(this%cell%atomic%moment/=0.0d0.and.this%cell%atomic%orbitals>0)) this%dim_modes(6)=this%dim_modes(6)*2 !another factor for spin
    endif
!    this%dim_modes(6)=this%dim_modes(6)*2   !UNIFY WITH HAMILTONIAN(MAGNETIC ORDER)
    !if(extpar_io%enable_u) dim_modes(5)=size(cell%atomic)*3
    this%order_set=this%dim_modes>0
    if(this%order_set(1)) this%site_per_cell(1)=this%nmag
    if(this%order_set(2)) this%site_per_cell(2)=1
    if(this%order_set(3)) this%site_per_cell(3)=1
    if(this%order_set(4)) this%site_per_cell(4)=1
    if(this%order_set(5)) this%site_per_cell(5)=this%nph
    if(this%order_set(6)) this%site_per_cell(6)=this%norb

    !new orderparameter format
    if(this%order_set(1)) Call this%M%init(this%dim_lat,this%dim_modes(1),dim_modes_inner(1),is_cmplx=.false.)
    if(this%order_set(2)) Call this%E%init(this%dim_lat,this%dim_modes(2),dim_modes_inner(2),is_cmplx=.false.)
    if(this%order_set(3)) Call this%B%init(this%dim_lat,this%dim_modes(3),dim_modes_inner(3),is_cmplx=.false.)
    if(this%order_set(4)) Call this%T%init(this%dim_lat,this%dim_modes(4),dim_modes_inner(4),is_cmplx=.false.)
    if(this%order_set(5)) Call this%u%init(this%dim_lat,this%dim_modes(5),dim_modes_inner(5),is_cmplx=.false.)
    if(this%order_set(6)) Call this%w%init(this%dim_lat,this%dim_modes(6),dim_modes_inner(6),is_cmplx=.true. )
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
    if(this%order_set(6).and..not.isinit(6)) Call this%w%read_file(trim(order_parameter_name(6))//suffix,isinit(6))
    deallocate(suffix)
    if(present(isinit_opt)) isinit_opt=isinit
end subroutine

subroutine lattice_position(this,ind,pos)
    !slow routine to get position of individual index
    class(lattice)          :: this
    integer,intent(in)      :: ind(4)  !(i_x,i_y,i_y,i_at)
    real(8),intent(out)     :: pos(3)

    real(8)                 :: ind_real(3)

    ind_real=real(ind(1:3)-1,8)
    pos=matmul(ind_real,this%areal)
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
    class(lattice),intent(in)               ::  this
    integer,intent(in)                      ::  order
    real(8),pointer,contiguous,intent(out)  ::  point(:)

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
    case(6)
        point=>this%w%all_modes
    case default
        write(error_unit,*) 'order:',order
        STOP 'failed to associate pointer in set_order_point'
    end select
end subroutine

subroutine set_order_point_inner(this,order,point)
    class(lattice),intent(in)               ::  this
    integer,intent(in)                      ::  order
    real(8),pointer,contiguous,intent(out)  ::  point(:,:)

    select case(order)
    case(1)
        point=>this%M%modes_in
    case(2)
        point=>this%E%modes_in
    case(3)
        point=>this%B%modes_in
    case(4)
        point=>this%T%modes_in
    case(5)
        point=>this%u%modes_in
    case(6)
        point=>this%w%modes_in
    case default
        write(error_unit,*) 'order:',order
        STOP 'failed to associate pointer in set_order_point'
    end select
end subroutine

function get_order_dim(this,order,ignore) result(dim_mode)
    class(lattice),intent(in)   :: this
    integer,intent(in)          :: order
    logical,intent(in),optional :: ignore
    integer                     :: dim_mode

#ifdef CPP_DEBUG
    !this could spot a programming mistake
    if(order<1.or.order>number_different_order_parameters)then
        write(error_unit,'(//A)') "Error trying to get order-dimension of unimplemented integer identifier."
        write(error_unit,'(A,I6)') "Order identifier:",order
        write(error_unit,'(A,I6)') "max. identifier :",number_different_order_parameters
        write(error_unit,'(A)') "This error is most probably caused by a badly implemented Hamiltonian"
        Error STOP "Check code"
    endif
#endif
    dim_mode=this%dim_modes(order)
    if(dim_mode<1)then
        if(present(ignore))then
            if(ignore) return
        endif
        write(error_unit,'(//2A)') "Error trying to access dimension of order parameter: ",trim(order_parameter_name(order))
        write(error_unit,'(A,I6)') "Its current order dimension is:", dim_mode
        if(dim_mode==0) write(error_unit,'(A,I6)') "Check input if the wanted order-parameter has been enabled"
        ERROR STOP "Check input"
    endif
end function

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
    if(this%order_set(6)) Call this%w%copy(copy%w,this%dim_lat)
end subroutine

subroutine copy_val_lattice(this,copy)
    class(lattice),intent(inout)    :: this,copy
   
    if(this%order_set(1)) Call this%M%copy_val(copy%M)
    if(this%order_set(2)) Call this%E%copy_val(copy%E)
    if(this%order_set(3)) Call this%B%copy_val(copy%B)
    if(this%order_set(4)) Call this%T%copy_val(copy%T)
    if(this%order_set(5)) Call this%u%copy_val(copy%u)
    if(this%order_set(6)) Call this%w%copy_val(copy%u)
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
    write(error_unit,'(//2A)') 'inserted name:',name_in
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
        write(error_unit,'(//A,I6)') "failing to get operator integer for ",abbrev_in(i:i)
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
    Call MPI_Bcast(this%areal        , 9                                , MPI_REAL8  , comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%astar        , 9                                , MPI_REAL8  , comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%dim_lat      , 3                                , MPI_INTEGER, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%Ncell        , 1                                , MPI_INTEGER, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%a_sc         , 9                                , MPI_REAL8  , comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%a_sc_inv     , 9                                , MPI_REAL8  , comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%dim_modes    , size(this%dim_modes)             , MPI_INTEGER, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%n_system     , 1                                , MPI_INTEGER, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%nmag         , 1                                , MPI_INTEGER, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%nph          , 1                                , MPI_INTEGER, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%periodic     , 3                                , MPI_LOGICAL, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%order_set    , size(this%order_set)             , MPI_LOGICAL, comm%mas, comm%com,ierr)
    Call MPI_Bcast(this%site_per_cell, number_different_order_parameters, MPI_INTEGER, comm%mas, comm%com,ierr)

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
    if(this%order_set(6)) Call this%w%bcast(comm)
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
    if(this%order_set(6)) Call this%w%bcast_val(comm)
#else
    continue
#endif
end subroutine

end module
