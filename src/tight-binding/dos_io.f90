module m_dos_io
use m_derived_types, only: lattice
use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
private
public dos_bnd_io , dos_orb_io
public dos_get_ind, dos_get_orb

type dos_bnd_io
    integer ::  atid=-1     !atom id
    integer ::  orb=-1      !orbital (0,atom%norb)
    integer ::  spin=-1     !spin (0,nspin)
    integer ::  site(3)=1   !site
contains
    procedure   :: dos_read => dos_read_bnd
    generic :: read(formatted) => dos_read
    procedure :: dos_assign => dos_assign_bnd
    generic :: assignment(=) => dos_assign
    procedure :: print_std => dos_print_bnd
    procedure :: check => dos_check_bnd
    procedure :: bcast => bcast_bnd
end type

type dos_orb_io
    integer ::  atid=-1     !atom id
    integer ::  orb=-1      !orbital (0,atom%norb)
    integer ::  spin=-1     !spin (1,2)
contains
    procedure   :: dos_read => dos_read_orb
    generic :: read(formatted) => dos_read
    procedure :: dos_assign => dos_assign_orb
    generic :: assignment(=) => dos_assign
    procedure :: print_std => dos_print_orb
    procedure :: check => dos_check_orb
    procedure :: bcast => bcast_orb
end type

contains


subroutine dos_get_orb(io,lat,nspin,norb_at,norb_at_off,orb_out)
    type(dos_orb_io),intent(in)         ::  io(:)
    type(lattice),intent(in)            ::  lat
    integer,intent(in)                  ::  nspin
    integer,intent(in)                  ::  norb_at(:)
    integer,intent(in)                  ::  norb_at_off(:)
    integer,intent(inout),allocatable   ::  orb_out(:)

    integer     :: i_cell, i_bnd
    integer     :: orb,spn
    integer     :: ndim

    ndim=sum(norb_at)*nspin
    allocate(orb_out(size(io)),source=0)
    do i_bnd=1,size(io)
        orb=io(i_bnd)%orb+norb_at_off(io(i_bnd)%atid)
        orb_out(i_bnd)=(orb-1)*nspin+io(i_bnd)%spin
    enddo
end subroutine

subroutine dos_get_ind(dos_bnd,lat,nspin,norb_at,norb_at_off,bnd_out)
    !subroutine which translates the dos_bnd_io to the start and end position of the contiguous
    !section in space of the Hamiltonian basis on which the projection is considered
    type(dos_bnd_io),intent(in)         ::  dos_bnd(:)
    type(lattice),intent(in)            ::  lat
    integer,intent(in)                  ::  nspin
    integer,intent(in)                  ::  norb_at(:)
    integer,intent(in)                  ::  norb_at_off(:)
    integer,intent(inout),allocatable   ::  bnd_out(:,:)

    integer     :: i_cell, i_bnd
    integer     :: orb(2),spn(2)
    integer     :: ndim

    ndim=sum(norb_at)*nspin

    allocate(bnd_out(2,size(dos_bnd)),source=0)
    do i_bnd=1,size(dos_bnd)
        i_cell=lat%index_m_1(dos_bnd(i_bnd)%site)
        orb=dos_bnd(i_bnd)%orb
        if(dos_bnd(i_bnd)%orb==0) orb=[1,norb_at(dos_bnd(i_bnd)%atid)]
        orb=orb+norb_at_off(dos_bnd(i_bnd)%atid)
        spn=dos_bnd(i_bnd)%spin
        if(dos_bnd(i_bnd)%spin==0) spn=[1,nspin]
        bnd_out(:,i_bnd)=(i_cell-1)*ndim+(orb-1)*nspin+spn
    enddo
end subroutine

subroutine dos_read_bnd(par, unit, iotype, v_list, iostat, iomsg)
    class(dos_bnd_io),intent(inout) :: par
    integer, intent(in)             :: unit
    character(*), intent(in)        :: iotype
    integer, intent(in)             :: v_list(:)
    integer, intent(out)            :: iostat
    character(*), intent(inout)     :: iomsg
    
    type(dos_bnd_io)                :: tmp
   
    read(unit,*,iostat=iostat,iomsg=iomsg) tmp%atid,tmp%orb,tmp%spin,tmp%site
    backspace(unit)
    if (iostat > 0)then    !try to read without z-site
        read(unit,*,iostat=iostat,iomsg=iomsg) tmp%atid,tmp%orb,tmp%spin,tmp%site(1:2)
        tmp%site(3)=1
        backspace(unit)
    endif
    if (iostat > 0)then    !try to read without z-site
        read(unit,*,iostat=iostat,iomsg=iomsg) tmp%atid,tmp%orb,tmp%spin,tmp%site(1)
        tmp%site(3)=1
        backspace(unit)
    endif
    read(unit,*)    !together with backspace make sure that is advances
    if(iostat==0) par=tmp
end subroutine

subroutine dos_read_orb(par, unit, iotype, v_list, iostat, iomsg)
    class(dos_orb_io),intent(inout) :: par
    integer, intent(in)             :: unit
    character(*), intent(in)        :: iotype
    integer, intent(in)             :: v_list(:)
    integer, intent(out)            :: iostat
    character(*), intent(inout)     :: iomsg
    
    type(dos_orb_io)                :: tmp
   
    read(unit,*,iostat=iostat,iomsg=iomsg) tmp%atid,tmp%orb,tmp%spin
    backspace(unit)
    read(unit,*)    !together with backspace make sure that is advances
    if(iostat==0) par=tmp
end subroutine


subroutine dos_assign_bnd(par,par_in)
    class(dos_bnd_io), intent(out):: par
    type(dos_bnd_io) , intent(in ):: par_in

    par%atid=par_in%atid  
    par%orb =par_in%orb
    par%spin=par_in%spin    
    par%site=par_in%site    
end subroutine

subroutine dos_assign_orb(par,par_in)
    class(dos_orb_io), intent(out):: par
    type(dos_orb_io) , intent(in ):: par_in

    par%atid=par_in%atid  
    par%orb =par_in%orb
    par%spin=par_in%spin    
end subroutine


subroutine dos_print_bnd(this,io_in)
    use, intrinsic :: iso_fortran_env, only : output_unit
    class(dos_bnd_io),intent(in)    :: this
    integer,intent(in),optional     :: io_in
    integer                         :: io_unit

    io_unit=output_unit
    if(present(io_in)) io_unit=io_in

    write(io_unit,'(A)')         'Local dos site-dependent input data:'
    write(io_unit,'(A, I6)')     '  atom id :', this%atid
    write(io_unit,'(A, I6)')     '  orbital :', this%orb
    write(io_unit,'(A, I6)')     '  spin    :', this%spin
    write(io_unit,'(A,3I6)')     '  site    :', this%site
end subroutine

subroutine dos_print_orb(this,io_in)
    use, intrinsic :: iso_fortran_env, only : output_unit
    class(dos_orb_io),intent(in)    :: this
    integer,intent(in),optional     :: io_in
    integer                         :: io_unit

    io_unit=output_unit
    if(present(io_in)) io_unit=io_in

    write(io_unit,'(A)')         'Local dos orbital-dependent input data:'
    write(io_unit,'(A, I6)')     '  atom id :', this%atid
    write(io_unit,'(A, I6)')     '  orbital :', this%orb
    write(io_unit,'(A, I6)')     '  spin    :', this%spin
end subroutine


subroutine dos_check_bnd(this,lat)
    class(dos_bnd_io), intent(in)   :: this
    type(lattice),intent(in)        :: lat

    if(this%atid<1.or.this%atid>size(lat%cell%atomic))then
        write(error_unit,'(//A/A)') "Atom id in local dos input out of range","Error in entry:" 
        Call this%print_std(error_unit)
        Stop "Check input"
    endif

    if(this%orb<0.or.this%orb>lat%cell%atomic(this%atid)%orbitals.or.lat%cell%atomic(this%atid)%orbitals==0)then
        write(error_unit,'(//A/A)') "Orbital in local dos out of range","Error in entry:" 
        Call this%print_std(error_unit)
        Stop "Check input"
    endif

    if(this%spin<0.or.this%spin>2)then 
        write(error_unit,'(//A/A)') "Spin in local dos out of range","Error in entry:" 
        Call this%print_std(error_unit)
        Stop "Check input"
    endif

    if(any(this%site<1).or.any(this%site>lat%dim_lat))then 
        write(error_unit,'(//A/A)') "site in local dos out of range","Error in entry:" 
        Call this%print_std(error_unit)
        Stop "Check input"
    endif
end subroutine

subroutine dos_check_orb(this,lat)
    class(dos_orb_io), intent(in)   :: this
    type(lattice),intent(in)        :: lat

    if(this%atid<1.or.this%atid>size(lat%cell%atomic))then
        write(error_unit,'(//A/A)') "Atom id in local dos input out of range","Error in entry:" 
        Call this%print_std(error_unit)
        Stop "Check input"
    endif

    if(this%orb<0.or.this%orb>lat%cell%atomic(this%atid)%orbitals.or.lat%cell%atomic(this%atid)%orbitals==0)then
        write(error_unit,'(//A/A)') "Orbital in local dos out of range","Error in entry:" 
        Call this%print_std(error_unit)
        Stop "Check input"
    endif

    if(this%spin<0.or.this%spin>2)then 
        write(error_unit,'(//A/A)') "Spin in local dos out of range","Error in entry:" 
        Call this%print_std(error_unit)
        Stop "Check input"
    endif
end subroutine

subroutine bcast_bnd(this,comm)
    use mpi_basic
    class(dos_bnd_io), intent(inout)    :: this
    type(mpi_type),intent(in)           :: comm
#ifdef CPP_MPI
    integer     :: val_arr(6), ierr
    val_arr(1)  =this%atid
    val_arr(2)  =this%orb
    val_arr(3)  =this%spin
    val_arr(4:6)=this%site
    Call MPI_bcast(val_arr, 6,MPI_INTEGER,comm%mas,comm%com,ierr)
    this%atid=val_arr(1)  
    this%orb =val_arr(2)  
    this%spin=val_arr(3)  
    this%site=val_arr(4:6)
#else
    continue
#endif
end subroutine

subroutine bcast_orb(this,comm)
    use mpi_basic
    class(dos_orb_io), intent(inout)    :: this
    type(mpi_type),intent(in)           :: comm
#ifdef CPP_MPI
    integer     :: val_arr(3), ierr
    val_arr(1)  =this%atid
    val_arr(2)  =this%orb
    val_arr(3)  =this%spin
    Call MPI_bcast(val_arr, 3,MPI_INTEGER,comm%mas,comm%com,ierr)
    this%atid=val_arr(1)  
    this%orb =val_arr(2)  
    this%spin=val_arr(3)  
#else
    continue
#endif
end subroutine

end module
