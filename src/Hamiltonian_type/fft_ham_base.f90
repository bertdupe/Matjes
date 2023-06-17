module m_fft_H_base
!module which contains the basic Hamiltonian description of a Hamiltonian solved through discrete fourier transformation
!while this type is not defined as abstract, in general it should be as its main set and evaluation routines will create a crash

!only the extensions of this type should be allocated

use m_type_lattice,only: lattice, dim_modes_inner 
use m_H_type, only: len_desc
use m_fft_H_internal, only: int_set_M, int_get_H
use, intrinsic :: iso_c_binding, only: C_int, C_DOUBLE
private
public fft_H

type            ::  fft_H
    private
    logical,public          :: set=.false.  !all parameters have been initialized
    integer(C_int),public   :: order=0      !lattice order parameter index (1=M,2=E etc.)
    integer(C_int),public   :: N_rep(3)=0   !size in each dimension of each considered field
    integer(C_int),public   :: dim_mode=0   !dimension of considered operator field
    logical,public          :: periodic(3)=.true.   !consider as periodic or open boundary condition along each direction (T:period, F:open)
                                            ! (dim_lat(i)=1->period(i)=T, since the calculation in the periodic case is easier, but choice of supercell_vec still does not consider periodicity)
    character(len=len_desc)     :: desc=""  !description of the Hamiltonian term, only used for user information and should be set manually 


    real(C_DOUBLE),allocatable,public  ::  M_n(:,:)    !magnetization in normal-space, set by M_internal
    real(C_DOUBLE),allocatable,public  ::  H_n(:,:)    !effective field normal-space, extract by H_internal

    procedure(int_set_M), pointer,nopass,public    ::  M_internal => null()    !function to set internal magnetization depending on periodic/open boundaries
    procedure(int_get_H), pointer,nopass,public    ::  H_internal => null()    !function to get effective field depending on periodic/open boundaries


contains
    !evaluation routines
    procedure,public    :: get_H            !get effective field
    procedure,public    :: get_H_single     !get effective field for a single site
    procedure,public    :: get_E            !get energy
    procedure,public    :: get_E_single     !get energy for a single site
    procedure,public    :: get_E_distrib    !get energy-distribution in Nmag*Ncell-space

    !initialization routines
    procedure,public    :: init_shape       !initializes the shape and the H and M operators, arrays
    procedure,public    :: init_op          !initializes the K_F array

    procedure,public    :: is_set           !returns set
    procedure,public    :: get_desc         !get the description
    
    !utility functions
    procedure,public    :: destroy          !destroys all data of this type
    procedure,public    :: copy             !copy this to new instance
    procedure,public    :: mv               !mv this to new instance
    procedure,public    :: add              !adds two fft_H in same space by adding their K_F
    procedure,public    :: bcast=>bcast_fft !mv this to new instance

    !internal procedures
    procedure           :: set_M                !set internal magnetization in normal-space from lattice
    procedure           :: init_internal        !initialize internal procedures
    procedure           :: same_space           !check if 2 fft_H act on same space

end type
contains 

subroutine get_desc(this,desc)
    class(fft_H),intent(in)             :: this
    character(len=len_desc),intent(out) :: desc

    desc=this%desc
end subroutine

function same_space(this,comp)result(same)
    class(fft_H),intent(in) :: this
    type(fft_H),intent(in)  :: comp
    logical                 :: same

    if(.not.this%set) ERROR STOP "CANNOT CHECK IF fft_H is the same space as this is not set"
    if(.not.comp%set) ERROR STOP "CANNOT CHECK IF fft_H is the same space as comp is not set"

    same=all(this%N_rep==comp%N_rep)
    same=same.and.all(this%periodic.eqv.comp%periodic)
    same=same.and.this%dim_mode==comp%dim_mode
end function

subroutine add(this,H_in)
    class(fft_H),intent(inout)  :: this
    class(fft_H),intent(in)     :: H_in

    if(.not.H_in%set) ERROR STOP "CANNOT ADD fft_H as H_in is not set"
    if(.not.this%set)then
        Call H_in%copy(this)
    else
        if(.not.this%same_space(H_in)) ERROR STOP "CANNOT ADD fft_H as this and H_in do not act on same space"
        this%desc="combined Hamiltonian"
    endif
end subroutine

subroutine bcast_fft(this,comm)
    use mpi_util
    class(fft_H),intent(inout)  ::  this
    type(mpi_type),intent(in)   ::  comm
    integer ::  shp(2)

    Call bcast(this%N_rep,comm)
    Call bcast(this%order,comm)
    Call bcast(this%periodic,comm)
    Call bcast(this%dim_mode,comm)
    Call bcast(this%desc,comm)

    if(comm%ismas) shp=shape(this%M_n)
    Call bcast(shp,comm)
    if(.not.comm%ismas)then
        Call this%init_internal(this%periodic)
        allocate(this%M_n(shp(1),shp(2))) 
        allocate(this%H_n(shp(1),shp(2))) 
    endif
end subroutine

subroutine mv(this,H_out)
    class(fft_H),intent(inout)  :: this
    class(fft_H),intent(inout)  :: H_out

    if(.not.this%set)then
        ERROR STOP "CANNOT MV UNINITIALIZED FFT_H"
    endif

    if(H_out%set)then
        ERROR STOP "CANNOT MV TO ALREADY INITIALIZED FFT_H"
    endif

    H_out%N_rep=this%N_rep
    H_out%order=this%order
    H_out%periodic=this%periodic
    H_out%dim_mode=this%dim_mode
    H_out%desc=this%desc

    if(allocated(this%M_N)) call move_alloc(this%M_n,H_out%M_n)
    if(allocated(this%H_N)) call move_alloc(this%H_n,H_out%H_n)

    H_out%M_internal=>this%M_internal
    H_out%H_internal=>this%H_internal

    H_out%set=this%set
end subroutine

subroutine copy(this,H_out)
    class(fft_H),intent(in)     :: this
    class(fft_H),intent(inout)  :: H_out

    if(.not.this%set)then
        ERROR STOP "CANNOT COPY UNINITIALIZED FFT_H"
    endif
    H_out%N_rep=this%N_rep
    H_out%order=this%order
    H_out%periodic=this%periodic
    H_out%dim_mode=this%dim_mode
    H_out%desc=this%desc

    allocate(H_out%M_n,mold=this%M_n)
    allocate(H_out%H_n,mold=this%H_n)

    H_out%M_internal=>this%M_internal
    H_out%H_internal=>this%H_internal

    H_out%set=this%set
end subroutine

subroutine destroy(this)
    class(fft_H),intent(inout)     :: this

    this%N_rep=0
    this%order=0
    this%set=.false.
    this%periodic=.true.
    this%dim_mode=0
    this%desc=""

    if(allocated(this%M_N)) deallocate(this%M_N)
    if(allocated(this%H_N)) deallocate(this%H_N)

    nullify(this%M_internal, this%M_internal)
end subroutine

pure function is_set(this)result(set)
    class(fft_H),intent(in)       :: this
    logical     ::  set
    set=this%set
end function

subroutine init_op(this,dim_mode,K_n,desc_in)
    !subroutine which initializes the fourier-transformed operator of K, while deallocating K_N
    class(fft_H),intent(inout)              :: this
    integer,intent(in)                      :: dim_mode
    real(8),intent(inout),allocatable       :: K_n(:,:,:)
    character(len=*),intent(in),optional    :: desc_in

    if(present(desc_in))then
        if(len(desc_in)<=len_desc) this%desc=desc_in
    endif
end subroutine

subroutine init_shape(this,abbrev,dim_mode,periodic,dim_lat,Kbd,N_rep)
    !initializes the arrays with whose the work will be done, sets the shapes, and returns shape data necessary for construction of the operator tensor
!$  use omp_lib    
    use m_type_lattice,only: op_abbrev_to_int
    class(fft_H),intent(inout)  :: this
    character(1),intent(in)     :: abbrev       !Abbreviation for order parameter (M,E,U,...) according to order_parameter_abbrev 
    integer,intent(in)          :: dim_mode
    logical,intent(in)          :: periodic(3)  !T: periodic boundary, F: open boundary
    integer,intent(in)          :: dim_lat(3)
    integer,intent(out)         :: Kbd(2,3) 
    integer,intent(out)         :: N_rep(3)

    integer         :: i
    integer(C_int)  :: Nk_tot           !number of state considered in FT (product of N_rep)
    integer(4)      :: tmp_order(1)

    N_rep=dim_lat
    tmp_order = op_abbrev_to_int(abbrev)
    this%order = tmp_order(1)
    !set K-boundaries for periodic boundaries
    Kbd(1,:)=[0,0,0]
    Kbd(2,:)=dim_lat-1
    !set K-boundaries for open boundaries if open
    do i=1,3
        if(.not.periodic(i))then
            N_rep(i)=2*N_rep(i)
            Kbd(:,i)=[-dim_lat(i)+1,dim_lat(i)-1]
        endif
    enddo
    this%N_rep=N_rep
    this%periodic=periodic
    this%dim_mode=dim_mode

    Nk_tot=product(N_rep)
    allocate(this%M_n(dim_mode,Nk_tot),source=0.0d0)
    allocate(this%H_n(dim_mode,Nk_tot),source=0.0d0)

    Call this%init_internal(periodic)
end subroutine

subroutine init_internal(this,periodic)
    !initialize internal procedures M_internal and H_internal
    use m_fft_H_internal
    class(fft_H),intent(inout)    :: this
    logical,intent(in)            :: periodic(3)  !T: periodic boundary, F: open boundary

    if(all(periodic))then
        this%M_internal=>set_M_period_TTT
        this%H_internal=>set_H_period_TTT
    elseif(all(periodic(1:2)))then
        this%M_internal=>set_M_period_TTF
        this%H_internal=>set_H_period_TTF
    elseif(all(periodic(1:1)))then
        this%M_internal=>set_M_period_TFF
        this%H_internal=>set_H_period_TFF
    else
        this%M_internal=>set_M_period_FFF
        this%H_internal=>set_H_period_FFF
    endif
end subroutine

subroutine get_E_distrib(this,lat,order,Htmp,E)
    class(fft_H),intent(inout)    ::  this
    type(lattice),intent(in)      ::  lat
    integer,intent(in)            ::  order     !index of considered order parameter
    real(8),intent(inout)         ::  Htmp(:,:)
    real(8),intent(out)           ::  E(:)

    !ugly implementation so that this works with all order parameters
    !(this should be cleaned up)
    real(8),pointer,contiguous    ::  M(:)
    real(8),pointer,contiguous    ::  M_v(:,:)
    integer(2)                    ::  shp(2)    !Shape, so that Htmps 2nd dimension indexes all sites 

    if(order == this%order) then     
        Call lat%set_order_point(this%order,M)
        M_v(1:size(Htmp,1),1:size(Htmp,2)) => M

        Call this%get_H(lat,Htmp)
        Htmp=Htmp*M_v
        shp=[dim_modes_inner(order),lat%site_per_cell(order)*lat%Ncell]
        E=sum(reshape(Htmp,shp),1)*2.0d0 !not sure about *2.0d0
    else
        E=0.0d0
    endif
end subroutine

subroutine get_E(this,lat,Htmp,E)
    class(fft_H),intent(inout)    ::  this
    type(lattice),intent(in)      ::  lat
    real(8),intent(inout)         ::  Htmp(:,:)
    real(8),intent(out)           ::  E

    !ugly implementation so that this works with all order parameters
    !(this should be cleaned up)
    real(8),pointer,contiguous    ::  M(:)
    real(8),pointer,contiguous    ::  M_v(:,:)
    
    Call lat%set_order_point(this%order,M)
    M_v(1:size(Htmp,1),1:size(Htmp,2)) => M

    Call this%get_H(lat,Htmp)
    Htmp=Htmp*M_v
    E=sum(Htmp)
end subroutine

subroutine get_E_single(this,lat,site,E)
    !get energy by getting the effective field caused by all sites and multiply the considered site's H with the moment there
    !it might be faster to do the explicit folding in real space to get the effective field caused by the considered site and multiply then with the entire magnetization (might be worth testing)
    class(fft_H),intent(inout)  :: this
    type(lattice),intent(in)    :: lat
    integer,intent(in)          :: site
    real(8),intent(out)         :: E

    real(8)                     :: Htmp(3)

    !ugly implementation so that this works with all order parameters
    !(this should be cleaned up)
    real(8),pointer,contiguous    ::  M(:)
    real(8),pointer,contiguous    ::  M_3(:,:)
    
    Call lat%set_order_point(this%order,M)
    M_3(1:3,1:size(M)/3) => M

    Call this%get_H_single(lat,site,Htmp)
    Htmp=Htmp*M_3(:,site)
    E=sum(Htmp)*2.0d0
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   Routines which have to be implemented locally    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine set_M(this,lat)
    !set this%M_n from lat%M%modes according to this%M_internal
    class(fft_H),intent(inout)    ::  this
    type(lattice),intent(in)      ::  lat

    ERROR STOP "set_M not implemented in base class"
end subroutine

subroutine get_H(this,lat,Hout)
    !get effective field H 
    class(fft_H),intent(inout)    ::  this
    type(lattice),intent(in)      ::  lat
    real(8),intent(inout)         ::  Hout(:,:)

    ERROR STOP "get_H not implemented in base class"
end subroutine

subroutine get_H_single(this,lat,site,Hout)
    !get effective field H
    class(fft_H),intent(inout)    ::  this
    type(lattice),intent(in)      ::  lat
    integer,intent(in)            ::  site
    real(8),intent(inout)         ::  Hout(3)

    ERROR STOP "get_H_single not implemented in base class"
end subroutine

end module
