module m_wannier_inp
    use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit, output_unit
private
public wann_dat

type wann_dat
    integer         :: num_wann=0
    integer         :: nrpts=0
    logical                 :: is_set=.false.
    logical                 :: spin_rearranged=.false.
    integer,allocatable     :: irvec(:,:)
    integer,allocatable     :: ndegen(:)
    complex(8),allocatable  :: H(:,:,:)
    !real(8),allocatable     :: r(:,:,:) !not necessary yet...
contains
    procedure :: read_file
    procedure :: combine_updn
    procedure :: rearrange_spin
    procedure :: destroy
end type

contains
subroutine destroy(this)
    class(wann_dat),intent(inout)  ::  this

    deallocate(this%irvec,this%ndegen,this%H)
    this%is_set=.false.
    this%spin_rearranged=.false.
    this%nrpts=0
    this%num_wann=0
end subroutine

subroutine combine_updn(this,up,dn)
    class(wann_dat),intent(inout)  ::  this
    type(wann_dat),intent(inout)   ::  up,dn
    integer         ::  num_wann_init

    if(.not.(up%is_set.and.dn%is_set)) return

    if(up%num_wann/=dn%num_wann)then
        write(error_unit,'(A,I6,A,I6)') "Cannot comine wannier up and down states as the number of wannier functions is inequivalent with:",up%num_wann," and",dn%num_wann
        STOP "ERROR"
    endif
    if(up%nrpts/=dn%nrpts)then
        write(error_unit,'(A,I6,A,I6)') "Cannot comine wannier up and down states as the number of R-vectors is inequivalent with:",up%nrpts," and",dn%nrpts
        STOP "ERROR"
    endif
    if(any(up%ndegen/=dn%ndegen))then
        write(error_unit,'(A)') "Cannot comine wannier up and down states as ndegen is inequivalent"
        STOP "ERROR"
    endif
    if(any(up%irvec/=dn%irvec))then
        write(error_unit,'(A)') "Cannot comine wannier up and down states as ndegen is inequivalent"
        STOP "ERROR"
    endif
    this%nrpts=up%nrpts
    this%num_wann=up%num_wann*2
    allocate(this%H(this%num_wann,this%num_wann,this%nrpts),source=(0.0d0,0.0d0))
    allocate(this%irvec,source=up%irvec)
    allocate(this%ndegen,source=up%ndegen)
    num_wann_init=up%num_wann
    this%H(1:num_wann_init,1:num_wann_init,:)=up%H
    this%H(num_wann_init+1:,num_wann_init+1:,:)=dn%H
    this%is_set=.true.
    Call up%destroy()
    Call dn%destroy()
    write(output_unit,'(A)') "Finished combining spin-up and spin-down Hamiltonians"
end subroutine

subroutine read_file(this,fname)
    !read wannier TB-input from wannier SEEDNAME_tb.dat file
    class(wann_dat),intent(inout)   ::  this
    character(*),intent(in)         ::  fname

    integer     ::  io
    real(8)     :: lat(3,3)
    integer     :: num_wann, nrpts
    integer     :: irpts
    logical     :: fexist
    integer     :: i,j

    character(10)   ::  char_nentry

    inquire(file=fname,exist=fexist)
    if(.not.fexist) return
    write(output_unit,'(/2A/)') "Reading wannier tight-binding input from: ",fname
    open(newunit=io,file=fname)
    read(io,*)
    read(io,*) lat
    lat=transpose(lat)
    read(io,*) num_wann
    read(io,*) nrpts
    this%num_wann=num_wann
    this%nrpts=nrpts
    allocate(this%ndegen(nrpts))
    read(io,'(15I5)') this%ndegen
    allocate(this%irvec(3,nrpts))
    allocate(this%H(num_wann,num_wann,nrpts),source=(0.0d0,0.0d0))
    write(char_nentry,'(I10)') num_wann*num_wann
    do irpts=1,nrpts !!
        read(io,'(/3I5)') this%irvec(:, irpts)
        do i=1,num_wann
            do j=1,num_wann
                read(io,'(13x,E15.8,1x,E15.8)') this%H(j,i,irpts)
            enddo
        enddo
    enddo
    !continue reading if real-space operator is necessary
    close(io)
    this%is_set=.true.
end subroutine

subroutine rearrange_spin(this)
    !wann input may be aranged with first spin-up and then spin-down states
    !this routine exchanges every second entry along the wannier directions to reorder
    !the Hamiltonian to the expected format where the spin if the innermost index
    class(wann_dat),intent(inout)   :: this
    complex(8),allocatable          :: tmp(:,:)
    integer                         :: shp(3)
    integer                         :: irpts,i
    integer                         :: wan_half


    if(this%spin_rearranged) STOP "Cannot rearrange spin of wannier, since this has already be done"
    this%spin_rearranged=.true.
    shp=shape(this%H)
    wan_half=shp(1)/2
    allocate(tmp(shp(1),shp(2)))
    do irpts=1,shp(3)
        tmp=this%H(:,:,irpts)
        do i=2,wan_half,2
            this%H(:,i,irpts)=tmp(:,i+wan_half)
            this%H(:,i+wan_half,irpts)=tmp(:,i)
        enddo
        tmp=transpose(this%H(:,:,irpts))
        do i=2,wan_half,2
            this%H(i,:,irpts)=tmp(:,i+wan_half)
            this%H(i+wan_half,:,irpts)=tmp(:,i)
        enddo
    enddo
end subroutine



end modulE

