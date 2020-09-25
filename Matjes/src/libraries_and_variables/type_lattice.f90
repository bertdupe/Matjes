module m_type_lattice
use m_basic_types

integer,parameter :: number_different_order_parameters=4    !m,E,B,T

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
     real(8) :: areal(3,3),astar(3,3),alat(3)
     integer :: dim_lat(3),n_system,dim_mode
     integer :: nmag !this nmag is nonsense and should be removed (number of m encoded elsewhere)
     integer, allocatable :: world(:)
     logical :: boundary(3)
! Table of pointer
     type(order_par)   :: ordpar !make this an array if using separated order parameters
     type(order_par)   :: M !magnetization 
     type(order_par)   :: E !Electric field
     type(order_par)   :: B !Magnetic field
     type(order_par)   :: T !Temperature
contains
    procedure :: copy => copy_lattice
    procedure :: copy_val_to => copy_val_lattice

end type lattice
private
public lattice, number_different_order_parameters

contains 

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
    Call self%ordpar%copy(copy%ordpar,self)
end subroutine


subroutine copy_val_lattice(self,copy)
    class(lattice),intent(inout)    :: self,copy
   
    !might make sense to check other values to make sure, but this should be faster
    Call self%ordpar%copy_val(copy%ordpar)
end subroutine

end module
