module m_delta_onsite
use m_types_tb_h_inp, only: TB_delta, TBio_delta_onsite_scf
use m_derived_types, only: lattice
private
public :: Hdelta

type Hdelta
    integer,allocatable     ::  orb(:,:)    
    complex(8),allocatable  ::  delta(:,:)
    real(8),allocatable     ::  V(:)
contains 

    procedure :: set
    procedure :: set_scf
    procedure :: print_file
    procedure :: read_file
    procedure :: bcast => bcast_local
end type

contains

subroutine bcast_local(this,comm)
    use mpi_basic
    use mpi_util
    class(Hdelta),intent(inout) :: this
    type(mpi_type),intent(in)   :: comm
#ifdef CPP_MPI
    Call bcast_alloc(this%orb  ,comm)
    Call bcast_alloc(this%delta,comm)
    Call bcast_alloc(this%V    ,comm)
#else
    continue
#endif
end subroutine 

subroutine set(this,delta,lat,norb_at_off)
    class(Hdelta),intent(inout)     :: this
    type(TB_delta),intent(in)       :: delta(:)       !delta parameter supplied input
    type(lattice),intent(in)        :: lat            !allknowing lattice input 
    integer,intent(in)              :: norb_at_off(:) !orbital offset at each atom

    integer                :: N_delta, i_delta
    integer                :: i, j
    integer,allocatable    :: ind_at(:) !atom indices for a given atom-type
    integer                :: mult

    N_delta=0
    do i=1,size(delta)
        if(delta(i)%dist/=0) STOP "delta distance /=0, so far only on-site delta terms implemented"
        if(delta(i)%attype(1)/=delta(i)%attype(2)) STOP "delta atomtypes differ, so far only on-site delta terms implemented"
        mult=1
!        if(delta(i)%orbital(1)/=delta(i)%orbital(2)) STOP "delta orbitals differ, so far only on-site delta terms implemented"
        if(delta(i)%orbital(1)/=delta(i)%orbital(2)) mult=2
        N_delta=N_delta+count(lat%cell%atomic(:)%type_id==delta(i)%attype(1))*mult
    enddo

    allocate(this%orb(2,N_delta),source=0)
    allocate(this%delta(lat%ncell,N_delta),source=(0.0d0,0.0d0))
    i_delta=0
    do i=1,size(delta)
        Call lat%cell%ind_attype(delta(i)%attype(1),ind_at)
        do j=1,size(ind_at)
            i_delta=i_delta+1
            this%orb(:,i_delta)=norb_at_off(ind_at(j))+delta(i)%orbital
            this%delta(:,i_delta)=delta(i)%val
            if(delta(i)%orbital(1)/=delta(i)%orbital(2))then
                i_delta=i_delta+1
                this%orb(:,i_delta)=norb_at_off(ind_at(j))+delta(i)%orbital(2:1:-1)
                this%delta(:,i_delta)=delta(i)%val
            endif
        enddo
        deallocate(ind_at)
    enddo
end subroutine


subroutine set_scf(this,scf_io,lat,norb_at_off)
    !set V-parameter to delta's already set
    !this is a rather cumbersome and error-prone way to set the self-consistent deltas with an initial guess, and it might be easier to get the entire delta and V from one input type
    class(Hdelta),intent(inout)             :: this
    type(TBio_delta_onsite_scf),intent(in)  :: scf_io(:)
    type(lattice),intent(in)                :: lat            !allknowing lattice input 
    integer,intent(in)                      :: norb_at_off(:) !orbital offset at each atom

    integer                 :: Nentry
    integer                 :: i,ii,iat
    integer                 :: nmult
    integer,allocatable     :: orb(:,:)
    real(8),allocatable     :: val(:)
    integer,allocatable     :: ind_at(:) !atom indices for a given atom-type

    if(.not.allocated(this%delta)) STOP "cannot set scf-delta if the deltas are not set before"

    !count how many scf entries there are
    allocate(this%V(size(this%orb,2)),source=0.0d0)
    Nentry=0
    do i=1,size(scf_io)
        nmult=1
        if(scf_io(i)%orbital(1)/=scf_io(i)%orbital(2)) nmult=2
        Nentry=Nentry+count(lat%cell%atomic(:)%type_id==scf_io(i)%attype)*nmult
    enddo

    !get all combinations of orbitals and attractive potential for onsite delta
    allocate(orb(2,Nentry),source=0)
    allocate(val(Nentry),source=0.0d0)
    ii=0
    do i=1,size(scf_io)
        Call lat%cell%ind_attype(scf_io(i)%attype,ind_at)
        do iat=1,size(ind_at)
            ii=ii+1
            orb(:,ii)=norb_at_off(ind_at(iat))+scf_io(i)%orbital
            val(ii)=scf_io(i)%val
            if(scf_io(i)%orbital(1)/=scf_io(i)%orbital(2))then
                ii=ii+1
                orb(:,ii)=norb_at_off(ind_at(iat))+scf_io(i)%orbital(2:1:-1)
                val(ii)=scf_io(i)%val
            endif
        enddo
        deallocate(ind_at)
    enddo

    !fill the set delta with an attractive potential, if supplied
    do i=1,size(orb,2)
        forall(ii = 1:size(this%orb,2), all(this%orb(:,ii)==orb(:,i))) this%V(ii)=val(i)
    enddo

    !do i=1,size(this%orb,2)
    !    write(*,*) i
    !    write(*,*) this%orb(:,i)
    !    write(*,*) this%delta(:,i)
    !    write(*,*) this%V(i)
    !    write(*,*)
    !enddo
end subroutine

subroutine print_file(this,fname,lat)
    class(Hdelta),intent(in)        :: this
    character(len=*),intent(in)     :: fname
    type(lattice),intent(in)        :: lat

    integer ::  io, i
    character(len=:),allocatable   :: frmt
    character(len=10)              :: flen

    open(newunit=io,file=fname)
    write(flen,'(I10)') size(this%delta,2)
    allocate(frmt,source='('//trim(adjustl(flen))//'(2E16.8,2X))')
    write(io,*) '!',lat%dim_lat, shape(this%delta)
    write(io,*) '!',this%orb
    do i=1,size(this%delta,1)
        write(io,frmt) this%delta(i,:)
    enddo
    close(io)
    deallocate(frmt)
end subroutine

subroutine read_file(this,fname,lat)
    use, intrinsic :: iso_fortran_env, only : error_unit
    class(Hdelta),intent(out)       :: this
    character(len=*),intent(in)     :: fname
    type(lattice),intent(in)        :: lat

    integer ::  io, i
    character(len=1)               :: tmp
    integer                        :: dim_lat(3)
    integer                        :: delta_shape(2)
    character(len=10)              :: flen
    character(len=:),allocatable   :: frmt

    open(newunit=io,file=fname)
    read(io,*) tmp,dim_lat, delta_shape
    if(any(lat%dim_lat/=dim_lat))then
        write(error_unit,'(2A)') "Cannot read onsite-delta from file: ", fname
        write(error_unit,'(A)') "dim_lat from file differs input system size"
        write(error_unit,'(A,3I6)') "system dim_lat=", lat%dim_lat
        write(error_unit,'(A,3I6)') "  read dim_lat=", dim_lat
        STOP
    endif
    allocate(this%orb(2,delta_shape(2)))
    read(io,*) tmp,this%orb
    if(any(this%orb>sum(lat%cell%atomic%orbitals)))then
        write(error_unit,'(2A)') "Cannot read onsite-delta from file: ", fname
        write(error_unit,'(A)')  "some read orbital is larger than than the number of orbitals in system"
        write(error_unit,'(A,I6)') "system sum orbitals=", sum(lat%cell%atomic%orbitals)
        write(error_unit,'(A,I6)') "  read max orbital =", maxval(this%orb)
        STOP
    endif
    allocate(this%delta(delta_shape(1),delta_shape(2)))
    write(flen,'(I10)') size(this%delta,2)
    allocate(frmt,source='('//trim(adjustl(flen))//'(2E16.8,2X))')
    do i=1,size(this%delta,1)
        read(io,frmt) this%delta(i,1)
    enddo
    close(io)
end subroutine



end module
