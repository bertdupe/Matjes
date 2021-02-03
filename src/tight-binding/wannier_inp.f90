module m_wannier_inp
private
public wann_dat

type wann_dat
    integer         :: num_wann=0
    integer         :: nrpts=0
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
    allocate(this%H(num_wann,num_wann,nrpts),source=(1.0d99,0.0d0))
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

end subroutine



end modulE

