module m_parameters_TB_IO_DOS
!module of type which contains the IO-parameters for the density of states (DOS) calculation
use m_dos_io
implicit none
private
public parameters_TB_IO_DOS

type parameters_TB_IO_DOS
    real(8)     :: E_ext(2)=[-1.0d0,1.0d0]          !minimum and maximum energy to plot in dos
    real(8)     :: dE=1.0d-2                        !energy binning size
    real(8)     :: sigma=1.0d-2                     !gauss smearing parameter for dos
    integer     :: kgrid(3)=[1,1,1]                 !number of k-points in each direction in case of k-space dos
    logical     :: print_kint=.false.               !print out the index of the currently considered k index 
    type(dos_bnd_io),allocatable    ::  bnd_io(:)   !io for local dos site dependent
    type(dos_orb_io),allocatable    ::  orb_io(:)   !io for local dos orbital dependent
    integer,allocatable :: bnd(:,:)                 !local dos bnd parameters (2,number local site dos)
    integer,allocatable :: orb(:)                   !local dos orbitals (number local orbital dos)
    character(:),allocatable    ::  fname_kmesh     !file name for kmesh integration grid input

    integer,allocatable :: fermi_orb(:)             !orbital indices of projections for fermi-surfaces
    logical             :: fermi_proj_all=.false.   !get fermi surface projection on all orbitals
contains
    procedure   :: read_file
    procedure   :: bcast => bcast_local
end type

contains

subroutine read_file(this,io,fname)
    use m_io_read_util
    use m_io_utils
    use, intrinsic :: iso_fortran_env, only : output_unit
    class(parameters_TB_IO_dos),intent(inout)   :: this
    integer,intent(in)                          :: io
    character(len=*), intent(in)                :: fname

    integer ::  N,ii,stat
    character(len=100) :: str

    call get_parameter(io,fname,'dos_sigma',this%sigma)
    call get_parameter(io,fname,'dos_E_ext',this%E_ext)
    call get_parameter(io,fname,'dos_dE',this%dE)
    call get_parameter(io,fname,'dos_kgrid',this%kgrid)
    call get_parameter(io,fname,'dos_print_kint',this%print_kint)

    str=" "
    Call get_parameter(io,fname,'dos_kmesh_file',str)
    if(len_trim(str)>1) this%fname_kmesh=trim(adjustl(str))

    N=0
    call get_parameter(io,fname,'TB_loc_dos',N)
    if(N>0)then
        Call set_pos_entry(io,fname,'TB_loc_dos')
        read(io,'(A)') str
        allocate(this%bnd_io(N))
        ii=1
        do while (ii<=size(this%bnd_io))
            read(io,'(a)',iostat=stat) str
            read(str,*,iostat=stat) this%bnd_io(ii)
            write(output_unit,'(A,I6,A)') 'TB_loc_dos entry no.',ii,':'
            Call this%bnd_io(ii)%print_std()
            ii=ii+1
        enddo 
        Call check_further_entry(io,fname,"TB_loc_dos")
    endif

    N=0
    call get_parameter(io,fname,'TB_orb_dos',N)
    if(N>0)then
        Call set_pos_entry(io,fname,'TB_orb_dos')
        read(io,'(A)') str
        allocate(this%orb_io(N))
        ii=1
        do while (ii<=size(this%orb_io))
            read(io,'(a)',iostat=stat) str
            read(str,*,iostat=stat) this%orb_io(ii)
            write(output_unit,'(A,I6,A)') 'TB_orb_dos entry no.',ii,':'
            Call this%orb_io(ii)%print_std()
            ii=ii+1
        enddo 
        Call check_further_entry(io,fname,"TB_orb_dos")
    endif

    N=0
    call get_parameter(io,fname,'N_fermi_orb',N)
    if(N/=0)then
        allocate(this%fermi_orb(N),source=0)
        call get_parameter(io,fname,'fermi_orb',N,this%fermi_orb)
    endif
    call get_parameter(io,fname,'fermi_proj_all',this%fermi_proj_all)
end subroutine

subroutine bcast_local(this,comm)
    use mpi_basic
    use mpi_util
    class(parameters_TB_IO_DOS),intent(inout)   ::  this
    type(mpi_type),intent(in)                   ::  comm
#ifdef CPP_MPI
    integer :: i,N
    logical :: used

    Call bcast(this%E_ext     ,comm)
    Call bcast(this%dE        ,comm)
    Call bcast(this%sigma     ,comm)
    Call bcast(this%kgrid     ,comm)
    Call bcast(this%print_kint,comm) 

    !site dependent dos input
    used=allocated(this%bnd_io)
    N=0
    if(used) N=size(this%bnd_io)
    Call bcast(N ,comm) 
    if(N>0)then
        if(.not.allocated(this%bnd_io)) allocate(this%bnd_io(N))
        do i=1,N
            Call this%bnd_io(i)%bcast(comm)
        enddo
    endif

    !orbital dependent dos input
    used=allocated(this%orb_io)
    if(used) N=size(this%orb_io)
    Call bcast(N ,comm) 
    if(N>0)then
        if(.not.allocated(this%orb_io)) allocate(this%orb_io(N))
        do i=1,N
            Call this%orb_io(i)%bcast(comm)
        enddo
    endif

    Call bcast_alloc(this%bnd,          comm)
    Call bcast_alloc(this%orb,          comm)
    Call bcast_alloc(this%fname_kmesh,  comm)
    Call bcast_alloc(this%fermi_orb,    comm)
    Call bcast(this%fermi_proj_all,     comm)

    Call bcast_alloc(this%fermi_orb, comm)
    Call bcast(this%fermi_proj_all, comm)
#else
    continue
#endif
end subroutine 

end module
