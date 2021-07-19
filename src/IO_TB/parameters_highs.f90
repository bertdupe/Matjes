module m_parameters_TB_IO_HIGHS
!module of type which contains the IO-parameters for the high-symmetry line calculation ( band structure)
implicit none
private
public parameters_TB_IO_HIGHS

type parameters_TB_IO_HIGHS
   integer                          ::  N_highsym=-1         !Number of high symmetry points to be connected
   real(8),allocatable              ::  k_highs(:,:)         !(3,Nighsym): high symmetry kpoints that are to be connected
   integer,allocatable              ::  k_highs_frac(:)      !(Nighsym): integer fractions to apply on input k_highs
   character(len=50),allocatable    ::  k_highs_label(:)     !high symmetry label names
   real(8)                          ::  aim_dist=1.0d-2      !aimed k-space distance between neighboring points along high symmetry line
contains
    procedure   :: read_file
    procedure   :: bcast => bcast_local
end type
contains

subroutine read_file(this,io,fname)
    use m_convert
    use m_io_read_util
    use m_io_utils
    class(parameters_TB_IO_highs),intent(inout) :: this
    integer,intent(in)                          :: io
    character(len=*), intent(in)                :: fname

    character(len=200)      :: string_tmp 
    character(len=50)       :: label_tmp
    integer                 :: i, stat
    logical                 :: success

    call get_parameter(io,fname,'N_highsym',this%N_highsym)
    if(this%N_highsym < 1)then
        return
    endif
    allocate(this%k_highs(3,this%N_highsym),source=0.0d0)
    allocate(this%k_highs_label(this%N_highsym))
    this%k_highs_label='"label"'
    allocate(this%k_highs_frac(this%N_highsym),source=1)
    call get_parameter(io,fname,'k_highs_frac',this%N_highsym,this%k_highs_frac)
    !read high symmetry positions
    Call set_pos_entry(io,fname,'k_highs_pts',success)
    read (io,*)
    do i=1,this%N_highsym
        read (io,'(a)',iostat=stat) string_tmp
        if (stat /= 0) ERROR STOP "Failed to read line after k_highs_pts, check input"
        read(string_tmp,*,iostat=stat) this%k_highs(:,i), this%k_highs_frac(i), label_tmp
        if(stat==0) this%k_highs_label(i)=label_tmp
        if(stat/=0)then
            read(string_tmp,*,iostat=stat) this%k_highs(:,i), this%k_highs_frac(i)
        endif
        if(stat/=0)then
            read(string_tmp,*,iostat=stat) this%k_highs(:,i)
        endif
        if(stat/=0) ERROR STOP "Failed to read k_highs_pts input, no 3 reals?"
    enddo
    call get_parameter(io,fname,'k_highs_dist',this%aim_dist)
end subroutine 

subroutine bcast_local(this,comm)
    use mpi_basic
    use mpi_util
    class(parameters_TB_IO_HIGHS),intent(inout) ::  this
    type(mpi_type),intent(in)                   ::  comm
#ifdef CPP_MPI
    Call bcast      (this%N_highsym           ,comm)
    Call bcast_alloc(this%k_highs             ,comm)
    Call bcast_alloc(this%k_highs_frac        ,comm)
    Call bcast_alloc(this%k_highs_label   ,50 ,comm)
    Call bcast      (this%aim_dist            ,comm)
#else
    continue
#endif
end subroutine 
end module
