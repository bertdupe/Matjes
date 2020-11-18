module m_io_gneb
use m_vector, only : norm
use m_type_lattice,only: lattice
implicit none

contains

subroutine read_inifin(file_ini,file_fin,images)
    use m_get_random
    use m_io_files_utils
    implicit none
    character(len=*), intent(in)    :: file_ini,file_fin
    type(lattice), intent(inout)    :: images(:)
    ! internal variables
     
    Call images(1)%M%read_file(file_ini)
    Call images(size(images))%M%read_file(file_fin)
end subroutine read_inifin


subroutine write_path(path)
use m_write_spin
use m_createspinfile
implicit none
real(kind=8), intent(in) :: path(:,:,:)
! internal variables
integer :: shape_path(3),i_nim

shape_path=shape(path)
do i_nim=1,shape_path(3)
   call WriteSpinAndCorrFile(i_nim,path(:,:,i_nim),'image-GNEB_')
   call CreateSpinFile('povray-GNEB_',i_nim,path(:,:,i_nim))
end do
end subroutine write_path


!
!
! Read the path from N files
!
!
subroutine read_path(fname_part,images,exists)
    use m_get_random
    use m_convert
    use m_io_files_utils
    implicit none
    character(len=*), intent(in) :: fname_part
    type(lattice), intent(inout)    :: images(:)
    logical, intent(out) :: exists
    ! internal variables
    logical ::  img_exist(size(images))
    integer :: i,i_nim,shape_path(3),io,N_cell
    real(kind=8) :: u(3),norm_local
    character(len=50) :: fname

    img_exist=.false.
    do i_nim = 1,size(images)
       fname=trim(fname_part)//'_'
       fname=convert(fname,i_nim)
       fname=trim(fname)//'.dat'
       io=open_file_read(fname)
       if (io.gt.0) then
          img_exist(i_nim)=.true.
          read(io,*) images(i)%M%modes_v
       else
          write(6,*) 'ERROR: File ',fname, ' does not exist. Path not loaded.'
       end if
       call close_file(fname,io)
    end do
    write(6,*) 'Path loaded.'
    exists=img_exist(1).and.img_exist(size(images)) !does exists only care if first & last file are read?

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
integer, intent(in) :: itr, itrmax
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
