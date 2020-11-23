module m_io_gneb
use m_vector, only : norm
use m_type_lattice,only: lattice
use m_input_types, only: GNEB_input
use m_convert, only: convert
use,intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT
implicit none
private
public :: read_path,read_path_inifin,prn_gneb_progress,write_en,write_path

contains

subroutine read_path_inifin(io_gneb,images)
    type(GNEB_input),intent(in)     :: io_gneb
    type(lattice), intent(inout)    :: images(:)

    logical     ::  exists
    character(:),allocatable        ::  momfile
    integer     :: nim

    inquire(file=io_gneb%momfile_i,exist=exists)
    if(exists)then
        write(OUTPUT_UNIT,'(2A)') "Reading initial GNEB image from: ",io_gneb%momfile_i
        Call images(1)%M%read_file(io_gneb%momfile_i)
    else
        momfile=convert(io_gneb%restartfile_if,'_1.dat')
        write(OUTPUT_UNIT,'(2A)') "Reading initial GNEB image from: ",momfile
        Call images(1)%M%read_file(momfile)
    endif
    if(allocated(momfile)) deallocate(momfile)

    nim=size(images)
    inquire(file=io_gneb%momfile_f,exist=exists)
    if(exists)then
        write(OUTPUT_UNIT,'(2A)') "Reading final GNEB image from: ",io_gneb%momfile_f
        Call images(nim)%M%read_file(io_gneb%momfile_f)
    else
        momfile=convert(io_gneb%restartfile_if,'_',nim,'.dat')
        write(OUTPUT_UNIT,'(2A)') "Reading final GNEB image from: ",momfile
        Call images(nim)%M%read_file(momfile)
    endif

end subroutine


subroutine write_path(images)
    use m_write_spin
    use m_createspinfile
    type(lattice), intent(in)    :: images(:)
    ! internal variables
    integer :: i_nim
    
    do i_nim=1,size(images)
       call WriteSpinAndCorrFile(i_nim,images(i_nim)%M%modes_v,'image-GNEB_')
       call CreateSpinFile('povray-GNEB_',i_nim,images(i_nim)%M%modes_v)
    end do
end subroutine write_path


!
!
! Read the path from N files
!
!
subroutine read_path(fname_part,images)
    use m_get_random
    use m_convert
    use m_io_files_utils
    implicit none
    character(len=*), intent(in) :: fname_part
    type(lattice), intent(inout)    :: images(:)
    ! internal variables
    logical ::  img_exist(size(images))
    integer :: i,i_im,shape_path(3),io,N_cell
    real(kind=8) :: u(3),norm_local
    character(len=50) :: fname

    img_exist=.false.
    do i_im = 1,size(images)
       fname=trim(fname_part)//'_'
       fname=convert(fname,i_im)
       fname=trim(fname)//'.dat'
       Call images(i_im)%M%read_file(fname)
    end do

end subroutine read_path



!> Print the path to file. Ensembles correspond to images in GNEB method
subroutine write_en(n,x,y,dy,x0,filn,do_norm_rx)
use m_io_files_utils
implicit none
integer, intent(in) :: n     !< Number of samples
real(kind=8), intent(in) :: x(n)        !< Reaction coordinate
real(kind=8), intent(in) :: y(n)        !< Energy
real(kind=8), intent(in) :: dy(n)       !< Derivative of the energy wrt x
real(kind=8), intent(in) :: x0        !< Normalization for the reaction coordinate
character(len=*), intent(in) :: filn             !< filename
character(len=1), intent(in) :: do_norm_rx   !< normalize reaction coordinate (Y/N)
! internal variables
integer :: i,io
real(kind=8) :: norm_local

if (do_norm_rx.eq.'Y') then
   norm_local = x0
elseif (do_norm_rx.eq.'N') then
   norm_local = 1d0
else
   write(6,*) 'Invalid value for do_norm_rx!'
   STOP
end if

io=open_file_write(filn)
do i=1,n
   write(io,'(3(E20.12E3,3x))') x(i)/norm_local, y(i), dy(i)
end do
call close_file(filn,io)

end subroutine write_en


subroutine prn_gneb_progress(itr,itrmax,fchk,imax,do_ci,ci)
implicit none
integer(8), intent(in) :: itr, itrmax
real(kind=8), intent(in) :: fchk
integer, intent(in) :: imax,ci
character(len=1), intent(in) :: do_ci

if (do_ci.eq.'Y') then
   write(6,'(a,2(E20.12E3,a),I8,a,I8)') 'MP  ',real(itr)/real(itrmax)*100d0,'% of itrmax.   fchk: ',fchk,'   imax: ',imax,'   ci: ',ci
else
   write(6,'(a,2(E20.12E3,a),I8)') 'MP  ',real(itr)/real(itrmax)*100d0,'% of itrmax.   fchk: ',fchk,'   imax: ',imax
end if

end subroutine prn_gneb_progress

end module m_io_gneb
