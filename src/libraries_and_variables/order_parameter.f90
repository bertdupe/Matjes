module m_order_parameter
implicit none
private
public order_par

type order_par
!variable that contains and order parameter field saved in data_real with modes/all_modes/modes_v pointing to that with different shapes

!data_real REQUIRE UTMOST CARE, DON'T CHANGE THEM EXTERNALLY AND THINK CAREFULLY ABOUT WHAT YOU DO
!UNDER NO CIRCUMSTANCES REMOVE PRIVATE FROM data_real
!THAT COULD LEAD TO A GIGANTIC MESS WITH MEMORY LEAKS

    real(8),pointer,private,contiguous      :: data_real(:)=>null() !main pointer that get allocated
    !pointers onto data_real, which can be adressed from outside
    real(8),pointer,contiguous              :: modes(:,:,:,:) => null()
    real(8),pointer,contiguous              :: all_modes(:)   => null()
    real(8),pointer,contiguous              :: modes_v(:,:)   => null() !(1:dimmode,:) shape
    real(8),pointer,contiguous              :: modes_3(:,:)   => null() !(1:3,:) shape, only associated if it makes sense (help accessing vector in 3-dimensional space)

    integer                                 :: dim_mode=0
contains 
    procedure :: init => init_order_par
    procedure :: copy => copy_order_par
    procedure :: copy_val => copy_val_order_par
    procedure :: delete => delete_order_par
    procedure :: read_file
    procedure :: truncate
    procedure :: bcast
    final :: final_order_par
end type
contains

subroutine bcast(this,comm)
use mpi_basic                
    class(order_par),intent(inout)  ::  this
    type(mpi_type),intent(in)       ::  comm

#ifdef CPP_MPI
    integer     :: ierr
    integer     :: N
    integer     :: mode_shape(4)
    integer     :: dim_lat(3)
    integer     :: dim_mode

    Call MPI_Bcast(this%dim_mode, 1, MPI_INTEGER, comm%mas, comm%com,ierr)
    if(comm%ismas)then
        N=size(this%data_real)
        mode_shape=shape(this%modes)
    endif
    Call MPI_Bcast(mode_shape, 4, MPI_INTEGER, comm%mas, comm%com,ierr)
    dim_mode=mode_shape(1)
    dim_lat=mode_shape(2:4)
    if(.not.comm%ismas)then
        Call this%init(dim_lat,dim_mode)
    endif
    Call MPI_Bcast(this%modes, size(this%modes), MPI_REAL8, comm%mas, comm%com,ierr)
#else
    continue
#endif
end subroutine

subroutine truncate(this,eps_in)
    class(order_par),intent(inout)      :: this
    real(8),intent(in),optional         :: eps_in
    real(8)     ::  eps
    real(8),parameter       ::  frac=1.0d-50

    if(present(eps_in))then
        eps=eps_in
    else
        eps=maxval(abs(this%all_modes))*frac
    endif
    where(abs(this%all_modes)<eps) this%all_modes=0.0d0
end subroutine

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

subroutine init_order_par(self,dim_lat,dim_mode,vec_val,val)
    !if the data pointers are not allocated, initialize them with 0
    !if the data pointers are allocated, check that size requirements are identical
    !associate public pointers to internal data storage
    class(order_par),intent(inout)  :: self
    integer,intent(in)              :: dim_lat(3)
    integer,intent(in)              :: dim_mode
    real(8),intent(in),optional     :: vec_val(dim_mode),val    !possibities to initialize order parameters
    !internal
    integer                         :: N,i
   
    !initialize real array data if necessary and associate pointers
    self%dim_mode=dim_mode
    N=dim_mode*product(dim_lat)
    if(.not.associated(self%data_real))then
        allocate(self%data_real(N),source=0.0d0)
    else
        if(size(self%data_real)/=N)then
            STOP "lattice shape does not match initializing already associated order_Par%modes"
        endif
    endif
    self%all_modes(1:N)=>self%data_real
    self%modes(1:dim_mode,1:dim_lat(1),1:dim_lat(2),1:dim_lat(3))=>self%data_real
    self%modes_v(1:dim_mode,1:product(dim_lat))=>self%data_real
    if(dim_mode>=3.and.modulo(dim_mode,3)==0) self%modes_3(1:3,1:product(dim_lat)*(dim_mode/3))=>self%data_real

    if(present(vec_val))then
        do i=1,size(self%modes_v,2)
            self%modes_v(:,i)=vec_val
        enddo
    endif

    if(present(val)) self%all_modes=val
end subroutine

subroutine final_order_par(self)
    type(order_par),intent(inout) :: self

    Call self%delete()
end subroutine

subroutine delete_order_par(self)
    class(order_par),intent(inout) :: self

    nullify(self%all_modes,self%modes,self%modes_v,self%modes_3)
    if(associated(self%data_real)) deallocate(self%data_real)
end subroutine

subroutine copy_order_par(self,copy,dim_lat)
    class(order_par),intent(in)    :: self
    class(order_par),intent(inout) :: copy
!    class(lattice),intent(in)      :: lat  !intent(in) might be a problem since a lattice might get destroyed <- check if finalization occurs
    integer,intent(in)              :: dim_lat(3)
   
    Call copy%delete()
    if(.not.associated(self%modes))then
        STOP "failed to copy order_par, since the source modes are not allocated"
    endif
    allocate(copy%data_real,source=self%data_real)
    Call copy%init(dim_lat,self%dim_mode)

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
