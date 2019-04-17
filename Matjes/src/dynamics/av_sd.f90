module m_sd_averages
use m_basic_types, only : vec_point
use m_topo_commons, only : Q_topo
interface sd_charge
   module procedure sd_charge_2D
!   module procedure sd_charge_2D_SL
end interface sd_charge
contains

function sd_charge_2D(i)
use m_vector, only : cross,norm,area
Implicit none
real(kind=8), dimension(4) :: sd_charge_2D
integer, intent(in) :: i
! internal
real(kind=8) :: vort(3),qtopo
integer :: N_spin,j

qtopo=0.0d0
vort=0.0d0

qtopo=area(Q_topo(1,i)%w,Q_topo(2,i)%w,Q_topo(3,i)%w)+area(Q_topo(1,i)%w,Q_topo(3,i)%w,Q_topo(4,i)%w)

vort=cross(Q_topo(1,i)%w,Q_topo(2,i)%w)+cross(Q_topo(2,i)%w,Q_topo(3,i)%w) &
        +cross(Q_topo(3,i)%w,Q_topo(4,i)%w)+cross(Q_topo(4,i)%w,Q_topo(1,i)%w)


sd_charge_2D(1)=qtopo
sd_charge_2D(2:4)=vort

end function sd_charge_2D

!      function sd_charge_2D_SL(i_x,i_y,i_z,i_m,spin,my_lattice)
!      use m_vector, only : cross,norm,area
!      Implicit none
!      type(lattice), intent(in) :: my_lattice
!      real(kind=8), dimension(4) :: sd_charge_2D_SL
!      integer, intent(in) :: i_x,i_y,i_z,i_m
!      real(kind=8), intent(in) :: spin(:,:,:,:,:)
!      integer :: ipu,ipv,ipm,dim_lat(3)
!      real(kind=8) :: vort(3),qtopo
!
!      dim_lat=my_lattice%dim_lat
!
!      ipu=mod(1+i_x+dim_lat(1)-1,dim_lat(1))+1
!      ipv=mod(i_y+1+dim_lat(2)-1,dim_lat(2))+1
!      ipm=mod(i_m+size(spin,5),size(spin,5))+1
!
!       qtopo=area(spin(:,i_x,i_y,i_z,i_m),spin(:,ipu,i_y,i_z,i_m),spin(:,ipu,ipv,i_z,i_m))+ &
!        area(spin(:,i_x,i_y,i_z,i_m),spin(:,ipu,ipv,i_z,i_m),spin(:,i_x,ipv,i_z,i_m))
!
!      vort=cross(Spin(:,i_x,i_y,i_z,i_m),Spin(:,ipu,i_y,i_z,i_m))+ &
!       cross(Spin(:,ipu,i_y,i_z,i_m),Spin(:,ipu,ipv,i_z,i_m)) &
!        +cross(Spin(:,ipu,ipv,i_z,i_m),Spin(:,i_x,ipv,i_z,i_m))+ &
!        cross(Spin(:,i_x,ipv,i_z,i_m),Spin(:,i_x,i_y,i_z,i_m))
!
!
!      sd_charge_2D_SL(1)=qtopo
!      sd_charge_2D_SL(2:4)=vort
!
!      end function sd_charge_2D_SL

end module m_sd_averages
