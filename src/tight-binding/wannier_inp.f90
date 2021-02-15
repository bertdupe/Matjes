module m_wannier_inp
private
public wann_dat

type wann_dat
    integer         :: num_wann=0
    integer         :: nrpts=0
    logical                 :: spin_rearranged=.false.
    integer,allocatable     :: irvec(:,:)
    integer,allocatable     :: ndegen(:)
    real(8),allocatable     :: r(:,:,:)
    complex(8),allocatable  :: H(:,:,:)
contains
    procedure :: read_file
end type

contains

subroutine read_file(this,fname)
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
    Call rearrange_spin(this)
end subroutine

subroutine rearrange_spin(this)
    !wann input may be aranged with first spin-up and then spin-down states
    !this routine exchanges every second entry along the wannier directions to reorder
    !the Hamiltonian to the expected format where the spin if the innermost index
    use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit
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
    write(error_unit,*) "CHECK rearrange_spin with bandstructure calculation"
end subroutine



end modulE

