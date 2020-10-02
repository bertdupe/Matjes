module m_type_lattice
use m_basic_types

integer,parameter :: number_different_order_parameters=4    !m,E,B,T

! unit cells
! the unit cell can be magnetic, ferroelectric or nothing and can also have transport
type t_cell
     type(atom), allocatable :: atomic(:)
     contains
     procedure :: ind_mag => cell_get_ind_mag
end type 

type order_par
!variable that contains all order parameters as vec_point and the
!corresponding target
!data_vecpoint and data_real REQUIRE UTMOST CARE, DON'T CHANGE THEM EXTERNALLY AND THINK WHAT YOU TO
!UNDER NO CIRCUMSTANCES REMOVE PRIVATE FROM data_real OR data_vecpoint
!THAT COULD LEAD TO A GIGANTIC MESS WITH MEMORY LEAK

    type(vec_point),pointer,private         :: data_vecpoint(:)=>null()
    real(8),pointer,private                 :: data_real(:)=>null()

    type(vec_point),pointer,contiguous      :: all_l_modes(:) => null()
    type(vec_point),pointer,contiguous      :: l_modes(:,:,:,:) => null() !main vec_type get allocated
    real(8),pointer,contiguous              :: modes(:,:,:,:,:) => null() !main pointer that get allocated
    real(8),pointer,contiguous              :: all_modes(:) => null()
    real(8),pointer,contiguous              :: modes_v(:,:) => null()

    integer                                 :: dim_mode=0
contains 
    procedure :: init => init_order_par
    procedure :: copy => copy_order_par
    procedure :: copy_val => copy_val_order_par
    procedure :: delete => delete_order_par
    final :: final_order_par
end type

! variable that defines the lattice
type lattice
     real(8)        :: areal(3,3),astar(3,3),alat(3)
     type(t_cell)   :: cell
     integer        :: dim_lat(3),n_system,dim_mode
     integer        :: nmag !this nmag is nonsense and should be removed (number of m encoded elsewhere)
     integer, allocatable :: world(:)
     logical :: boundary(3)

!order parameters
     type(order_par)   :: M !magnetization 
     type(order_par)   :: E !Electric field
     type(order_par)   :: B !Magnetic field
     type(order_par)   :: T !Temperature
     !ultimately delete ordpar
     type(order_par)   :: ordpar !make this an array if using separated order parameters

contains
    procedure :: copy => copy_lattice
    procedure :: get_pos_mag => lattice_get_position_mag
    procedure :: copy_val_to => copy_val_lattice
    procedure :: get_order_dim  => get_order_dim
    procedure :: set_order_point => set_order_point
    procedure :: index_m_1 => index_m_1
    procedure :: index_1_3 => index_1_3
    procedure :: used_order => lattice_used_order
end type lattice


private
public lattice, number_different_order_parameters, t_cell

contains 

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
    indm(3)=ind1/prod
    tmp=ind1-prod*indm(3)
    prod=this%dim_lat(1)
    indm(2)=tmp/prod
    indm(1)=tmp-indm(2)*prod
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
        STOP "dim_mode not positive, requested order parameter not set?"
    endif

end function

subroutine cell_get_ind_mag(this,ind_Nat)
    class(t_cell),intent(in)    ::  this
    integer                     ::  Nat
    integer                     ::  i,j
    integer,allocatable         ::  ind_Nat(:)

    Nat=0
    do i=1,size(this%atomic)
        if(this%atomic(i)%moment/=0.0d0) Nat=Nat+1
    enddo
    allocate(ind_Nat(Nat),source=0)
    j=0
    do i=1,size(this%atomic)
        if(this%atomic(i)%moment/=0.0d0)then
            j=j+1
            ind_Nat(j)=i
        endif
    enddo

end subroutine

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


subroutine init_order_par(self,lat,dim_mode)
    !if the data pointers are not allocated, initialize them with 0
    !if the data pointers are allocated, check that size requirements are identical
    !associate public pointers to internal data storage
    class(order_par),intent(inout) :: self
    class(lattice),intent(in)    :: lat
    integer,intent(in)           :: dim_mode
    integer                      :: N
    integer                      :: l,k,j,i          
    integer                      :: shape_modes(5),shape_lmodes(4)
    
    !initialize real array data if necessary and associate pointers
    self%dim_mode=dim_mode
    N=dim_mode*product(lat%dim_lat)*lat%nmag
    if(.not.associated(self%data_real))then
        allocate(self%data_real(N),source=0.0d0)
    else
        if(size(self%data_real)/=N)then
            STOP "lattice shape does not match initializing already associated order_Par%modes"
        endif
    endif
    self%all_modes(1:N)=>self%data_real
    self%modes(1:dim_mode,1:lat%dim_lat(1),1:lat%dim_lat(2),1:lat%dim_lat(3),1:lat%nmag)=>self%data_real
    self%modes_v(1:dim_mode,1:product(lat%dim_lat)*lat%nmag)=>self%data_real

    !initialize vec_pointer data if necessary and associate pointers
    N=product(lat%dim_lat)*lat%nmag
    if(.not.associated(self%data_vecpoint))then
        !allocate(self%data_vecpoint(lat%dim_lat(1),lat%dim_lat(2),lat%dim_lat(3),lat%nmag))
        allocate(self%data_vecpoint(N))
    else
        if(size(self%data_vecpoint)/=N)then
            STOP "lattice shape does not match initializing already allocated order_par%data_vecpoint"
        endif
    endif
    self%l_modes(1:lat%dim_lat(1),1:lat%dim_lat(2),1:lat%dim_lat(3),1:lat%nmag)=>self%data_vecpoint
    self%all_l_modes(1:N)=>self%data_vecpoint
    do l=1,lat%nmag
       do k=1,lat%dim_lat(3)
          do j=1,lat%dim_lat(2)
             do i=1,lat%dim_lat(1)
                self%l_modes(i,j,k,l)%w=>self%modes(:,i,j,k,l)
             enddo
          enddo
       enddo
    enddo
end subroutine

subroutine final_order_par(self)
    type(order_par),intent(inout) :: self

    Call self%delete()
end subroutine

subroutine delete_order_par(self)
    class(order_par),intent(inout) :: self
    
    nullify(self%all_l_modes,self%l_modes)
    if(associated(self%data_vecpoint)) deallocate(self%data_vecpoint)

    nullify(self%all_modes,self%modes)
    if(associated(self%data_real)) deallocate(self%data_real)
end subroutine

subroutine copy_order_par(self,copy,lat)
    class(order_par),intent(inout) :: self
    class(order_par),intent(inout) :: copy
    class(lattice),intent(in)      :: lat  !intent(in) might be a problem since a lattice might get destroyed <- check if finalization occurs
   
    Call copy%delete()
    if(.not.associated(self%modes))then
        STOP "failed to copy order_par, since the source modes are not allocated"
    endif
    allocate(copy%data_real,source=self%data_real)
    Call copy%init(lat,self%dim_mode)

end subroutine


subroutine copy_val_order_par(self,copy)
    class(order_par),intent(inout) :: self
    class(order_par),intent(inout) :: copy

   
    if(.not.associated(self%data_real))then
        STOP "failed to copy_val order_par, since the source modes are not allocated"
    endif
    if(.not.associated(copy%data_real))then
        STOP "failed to copy_val order_par, since the target modes are not allocated"
    endif

    if(size(self%data_real)/=size(copy%data_real))then
        STOP "failed to copy_val order_par, since their shapes are inconsistent"
    endif

    copy%all_modes=self%all_modes

end subroutine


subroutine copy_lattice(self,copy)
    class(lattice),intent(inout)    :: self,copy
    
    copy%areal=self%areal
    copy%astar=self%astar
    copy%alat=self%alat
    copy%dim_lat=self%dim_lat
    copy%n_system=self%n_system
    copy%dim_mode=self%dim_mode
    copy%nmag=self%nmag
    copy%boundary=self%boundary
    if(allocated(self%world)) allocate(copy%world,source=self%world)
    if(self%ordpar%dim_mode>0) Call self%ordpar%copy(copy%ordpar,self)

    if(self%m%dim_mode>0) Call self%M%copy(copy%M,self)
    if(self%e%dim_mode>0) Call self%E%copy(copy%E,self)
    if(self%b%dim_mode>0) Call self%B%copy(copy%B,self)
    if(self%t%dim_mode>0) Call self%T%copy(copy%T,self)
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

end module
