subroutine GNEB(my_lattice,my_motif,io_simu,ext_param)
use m_gneb_parameters
use m_gneb_utils
use m_write_spin
use m_createspinfile
use m_derived_types
use m_initialize_path
use m_torques, only : get_torques

implicit none
type(lattice), intent(in) :: my_lattice
type(cell), intent(in) :: my_motif
type(io_parameter), intent(in) :: io_simu
type(simulation_parameters), intent(in) :: ext_param
!internal variable
real(kind=8), allocatable :: path(:,:,:),spinsp(:,:)
real(kind=8), allocatable :: rx(:),ene(:),dene(:),xx(:),yy(:),dyy(:),c(:,:)
real(kind=8) :: rx0(1),drx0,ene0(1),dene0(1)
integer :: imax,i,ci,size_order,N_cell
logical :: gra_log

N_cell=product(my_lattice%dim_lat)
size_order=my_lattice%dim_mode
gra_log=io_simu%io_Xstruct
ci=1

call set_gneb_defaults()
call read_gneb_parameters()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize the different torques
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call get_torques('input')

allocate(path(size_order,N_cell,nim),spinsp(size_order,N_cell),rx(nim),ene(nim),dene(nim))
path=0.0d0
spinsp=0.0d0
rx=0.0d0
ene=0.0d0
dene=0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! initialization of the path
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call path_initialization(path,io_simu)

call write_path(path)
      
if (do_gneb=='Y') then
   write (6,'(a)') "GNEB calculation in progress..."
   call find_path(nim,N_cell,vpodt,vpomass,spring,mepftol,mepitrmax,meptraj_step,rx,ene,dene,path,my_lattice,io_simu)
         
   write (6,'(a)') "Done!"
end if
      
!if (do_gneb_ci=='Y') then
!   write (6,'(a)') "CI-GNEB calculation in progress..."
!   call find_path_ci(nim,N_cell,vpodt,vpomass,spring,path,mepftol_ci,mepitrmax,meptraj_step,my_lattice,rx,ene,dene,ci)
!   write(6,'(a)') *,'ci:',ci
!   write (6,'(a)') "Done!"
!end if
      
call write_en(nim,rx,ene,dene,rx(nim),'en_path.out',do_norm_rx)
call write_path(path)
      
allocate(xx(sample_num),yy(sample_num),dyy(sample_num),c(4,nim))
xx=0.0d0
yy=0.0d0
dyy=0.0d0
c=0.0d0

call hermite_fit(nim,sample_num,rx,ene,dene,xx,yy,dyy,c)

call write_en(sample_num,xx,yy,dyy,rx(nim),'enfit_path.out',do_norm_rx)
      
      
if (do_gneb_ci=='Y') then
   rx0(1) = rx(ci)
   do i=1,N_cell
      spinsp(:,i) = path(:,i,ci)
   end do
else
      
   call find_SP(nim,rx,c,imax,drx0)
   rx0(1) = rx(imax)+drx0
   if (imax==nim) then
      do i=1,N_cell
         spinsp(:,i) = path(:,i,imax)
      end do
   else
      call find_SP_conf(path(:,:,imax),path(:,:,imax+1),rx(imax),rx(imax+1),rx(imax)+drx0,spinsp)
   end if
end if

call spline_hermite_val(nim,rx,c,rx0(1),ene0(1),dene0(1))
      
call write_en(1,rx0,ene0,dene0,rx(nim),'ensp_path.out',do_norm_rx)
      
      
      
call WriteSpinAndCorrFile(spinsp,'image-GNEB-saddlepoint.dat')
call CreateSpinFile(spinsp,'povray-GNEB-saddlepoint.dat')
      

deallocate(path,rx,ene,dene,xx,yy,dyy,c)

end subroutine GNEB 
