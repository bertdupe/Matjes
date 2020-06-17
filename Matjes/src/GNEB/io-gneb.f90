module m_io_gneb
use m_vector, only : norm

contains

subroutine read_inifin(file_ini,file_fin,amp_rnd,path)
use m_get_random
use m_io_files_utils
implicit none
character(len=*), intent(in) :: file_ini,file_fin
real(kind=8), intent(in):: amp_rnd
real(kind=8), intent(inout) :: path(:,:,:)
! internal variables
integer :: i,io_ini,io_fin,N_cell,nim,shape_path(3)
real(kind=8) :: u(3),norm_local

shape_path=shape(path)
N_cell=shape_path(2)
nim=shape_path(3)

io_ini=open_file_read(file_ini)
if (io_ini.lt.0) then
    write(6,*) 'no input file for the initial state!'
    STOP
endif

io_fin=open_file_read(file_fin)
if (io_fin.lt.0) then
   write(6,*) 'no input file for the final state!'
   STOP
endif

do i=1,N_cell
   read(io_ini,*)  path(:,i,1)
end do
call close_file(file_ini,io_ini)

do i=1,N_cell
    read(io_fin,*) path(:,i,nim)
end do
call close_file(file_fin,io_fin)

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
! Read the path into N files
!
!
subroutine read_path(fname_part,amp_rnd,path,exists)
use m_get_random
use m_convert
use m_io_files_utils
implicit none
character(len=*), intent(in) :: fname_part
real(kind=8), intent(in) :: amp_rnd
real(kind=8), intent(inout) :: path(:,:,:)
logical, intent(out) :: exists
! internal variables
integer :: i,i_nim,shape_path(3),io,N_cell,nim
real(kind=8) :: u(3),norm_local
character(len=50) :: fname

shape_path=shape(path)
nim=shape_path(3)
N_cell=shape_path(2)

do i_nim = 1,nim
   fname=convert(fname_part,'_',i_nim,'.dat')

   io=open_file_read(fname)
   if (io.gt.0) then
      exists=.true.
      do i=1,N_cell
         read(io,*) path(:,i,i_nim)
      end do
   else
      exists=.false.
      write(6,*) 'ERROR: File ',fname, ' does not exist. Path not loaded.'
   end if
   call close_file(fname,io)
end do

write(6,*) 'Path loaded.'

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
