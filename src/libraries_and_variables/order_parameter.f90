module m_order_parameter
use m_basic_types, only: vec_point
implicit none
private
public order_par

! unit cells
! the unit cell can be magnetic, ferroelectric or nothing and can also have transport

type order_par
!variable that contains all order parameters as vec_point and the
!corresponding target
!data_vecpoint and data_real REQUIRE UTMOST CARE, DON'T CHANGE THEM EXTERNALLY AND THINK WHAT YOU TO
!UNDER NO CIRCUMSTANCES REMOVE PRIVATE FROM data_real OR data_vecpoint
!THAT COULD LEAD TO A GIGANTIC MESS WITH MEMORY LEAK

    type(vec_point),pointer,private,contiguous         :: data_vecpoint(:)=>null()
    real(8),pointer,private,contiguous                 :: data_real(:)=>null()

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
    procedure :: read_file
    final :: final_order_par
end type
contains

subroutine read_file(this,fname,fexist)
    use m_io_files_utils, only: open_file_read,close_file
    use,intrinsic ::  ISO_FORTRAN_ENV ,only: ERROR_UNIT,OUTPUT_UNIT
    class(order_par),intent(inout)      :: this
    character(*),intent(in)             :: fname
    logical,intent(out),optional        :: fexist

    integer     :: io,stat
    real(8)     :: tst

    if(present(fexist))then
        inquire(file=fname,exist=fexist)
        if(.not.fexist) return
    endif
    io=open_file_read(fname)
    if (io.lt.0) then
        write(ERROR_UNIT,'(//2A)') 'Failed to open file:',fname
        ERROR STOP "ERROR READING ORDER PARAMETER"
    endif
    read(io,*,iostat=stat)  this%modes_v
    if(stat>0)then
        write(ERROR_UNIT,'(//2A)') 'ERROR READING FILE:',fname
        write(ERROR_UNIT,'(2A)') 'Unexpected characters?'
        ERROR STOP "ERROR READING ORDER PARAMETER"
    elseif(stat<0)then
        write(ERROR_UNIT,'(//2A)') 'UNEXPECTED END OF FILE:',fname
        write(ERROR_UNIT,'(2A)') 'incompatible lattice size between input and file?'
        ERROR STOP "ERROR READING ORDER PARAMETER"
    else
        read(io,*,iostat=stat)  tst
        if(stat==0)then
            write(ERROR_UNIT,'(//2A)') 'UNEXPECTED FURTHER ENTRIES IN FILE:',fname
            write(ERROR_UNIT,'(2A)') 'incompatible lattice size between input and file?'
            ERROR STOP "ERROR READING ORDER PARAMETER"
        endif
    endif
    write(OUTPUT_UNIT,'(/2A/)') "Read input order parameter from: ",fname
    call close_file(fname,io)
end subroutine

subroutine init_order_par(self,dim_lat,nmag,dim_mode,vec_val,val)
    !if the data pointers are not allocated, initialize them with 0
    !if the data pointers are allocated, check that size requirements are identical
    !associate public pointers to internal data storage
    class(order_par),intent(inout)  :: self
    integer,intent(in)              :: dim_lat(3)
    integer,intent(in)              :: nmag
!    class(lattice),intent(in)    :: lat
    integer,intent(in)              :: dim_mode
    real(8),intent(in),optional     :: vec_val(dim_mode),val    !possibities to initialize order parameters
    integer                         :: N
    integer                         :: l,k,j,i          
   
    if(nmag>1) ERROR STOP "THIS STILL USES OLD SETIP WITH NMAG AS LAST PARAMETER"
    !initialize real array data if necessary and associate pointers
    self%dim_mode=dim_mode
    N=dim_mode*product(dim_lat)*nmag
    if(.not.associated(self%data_real))then
        allocate(self%data_real(N),source=0.0d0)
    else
        if(size(self%data_real)/=N)then
            STOP "lattice shape does not match initializing already associated order_Par%modes"
        endif
    endif
    self%all_modes(1:N)=>self%data_real
    self%modes(1:dim_mode,1:dim_lat(1),1:dim_lat(2),1:dim_lat(3),1:nmag)=>self%data_real
    self%modes_v(1:dim_mode,1:product(dim_lat)*nmag)=>self%data_real

    if(present(vec_val))then
        do i=1,size(self%modes_v,2)
            self%modes_v(:,i)=vec_val
        enddo
    endif

    if(present(val)) self%all_modes=val

    !initialize vec_pointer data if necessary and associate pointers
    N=product(dim_lat)*nmag
    if(.not.associated(self%data_vecpoint))then
        !allocate(self%data_vecpoint(dim_lat(1),dim_lat(2),dim_lat(3),nmag))
        allocate(self%data_vecpoint(N))
    else
        if(size(self%data_vecpoint)/=N)then
            STOP "lattice shape does not match initializing already allocated order_par%data_vecpoint"
        endif
    endif
    self%l_modes(1:dim_lat(1),1:dim_lat(2),1:dim_lat(3),1:nmag)=>self%data_vecpoint
    self%all_l_modes(1:N)=>self%data_vecpoint
    do l=1,nmag
       do k=1,dim_lat(3)
          do j=1,dim_lat(2)
             do i=1,dim_lat(1)
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

subroutine copy_order_par(self,copy,dim_lat,nmag)
    class(order_par),intent(in)    :: self
    class(order_par),intent(inout) :: copy
!    class(lattice),intent(in)      :: lat  !intent(in) might be a problem since a lattice might get destroyed <- check if finalization occurs
    integer,intent(in)              :: dim_lat(3),nmag
   
    Call copy%delete()
    if(.not.associated(self%modes))then
        STOP "failed to copy order_par, since the source modes are not allocated"
    endif
    allocate(copy%data_real,source=self%data_real)
    Call copy%init(dim_lat,nmag,self%dim_mode)

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
end module
