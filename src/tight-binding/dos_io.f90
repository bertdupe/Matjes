module m_dos_io
use m_derived_types, only: lattice
use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
private
public dos_bnd_io, dos_get_ind

type dos_bnd_io
    integer ::  atid=-1     !atom id
    integer ::  orb=-1      !orbital (0,atom%norb)
    integer ::  spin=-1     !spin (0,nspin)
    integer ::  site(3)=1   !site
contains
    procedure   :: dos_read
    generic :: read(formatted) => dos_read
    procedure :: dos_assign
    generic :: assignment(=) => dos_assign
    procedure :: print_std => dos_print
    procedure :: check => dos_check
end type
contains

subroutine dos_get_ind(dos_bnd,lat,nspin,norb_at,norb_at_off,bnd_out)
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

subroutine dos_read(par, unit, iotype, v_list, iostat, iomsg)
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

subroutine dos_assign(par,par_in)
    class(dos_bnd_io), intent(out):: par
    type(dos_bnd_io) , intent(in ):: par_in

    par%atid=par_in%atid  
    par%orb =par_in%orb
    par%spin=par_in%spin    
    par%site=par_in%site    
end subroutine

subroutine dos_print(this,io_in)
    use, intrinsic :: iso_fortran_env, only : output_unit
    class(dos_bnd_io),intent(in)    :: this
    integer,intent(in),optional     :: io_in
    integer                         :: io_unit

    io_unit=output_unit
    if(present(io_in)) io_unit=io_in

    write(io_unit,'(A)')         'Local dos input data:'
    write(io_unit,'(A, I6)')     '  atom id :', this%atid
    write(io_unit,'(A, I6)')     '  orbital :', this%orb
    write(io_unit,'(A, I6)')     '  spin    :', this%spin
    write(io_unit,'(A,3I6)')     '  site    :', this%site
end subroutine

subroutine dos_check(this,lat)
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

end module
