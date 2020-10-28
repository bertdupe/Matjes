module m_topo_commons
use m_derived_types, only : operator_real,t_cell,lattice
use m_basic_types, only : vec_point
Implicit none

interface get_charge
  module procedure :: get_charge_2D   ! take only one int argument (iomp)
  module procedure :: get_all_charge
  module procedure :: get_charge_lat
  module procedure :: get_charge_lat_i
end interface

type(vec_point),allocatable,protected,save,public :: Q_topo(:,:)

private
public :: get_size_Q_operator,associate_Q_operator,get_charge,neighbor_Q

contains

!
! function that returns the integral of the topological charge of the magnetic textures
!
function get_all_charge()
    Implicit none
    real(kind=8), dimension(5) :: get_all_charge
    ! internal
    integer :: i,N
    
    get_all_charge=0.0d0
    N=size(Q_topo,2)
    
    do i=1,N
       get_all_charge=get_all_charge+get_charge(i)
    enddo
    
end function get_all_charge

!
! function that returns the topological charge of the magnetic textures
!
function get_charge_2D(i)
use m_vector, only : cross,norm,area
Implicit none
real(kind=8), dimension(5) :: get_charge_2D
integer, intent(in) :: i
! internal
real(kind=8) :: vort(3),qtopo

qtopo=0.0d0
vort=0.0d0
get_charge_2D=0.0d0

qtopo=area(Q_topo(1,i)%w,Q_topo(2,i)%w,Q_topo(3,i)%w)+area(Q_topo(1,i)%w,Q_topo(3,i)%w,Q_topo(4,i)%w)

vort=cross(Q_topo(1,i)%w,Q_topo(2,i)%w,1,3)+cross(Q_topo(2,i)%w,Q_topo(3,i)%w,1,3) &
        +cross(Q_topo(3,i)%w,Q_topo(4,i)%w,1,3)+cross(Q_topo(4,i)%w,Q_topo(1,i)%w,1,3)

if (qtopo.gt.0.0d0) get_charge_2D(1)=qtopo
if (qtopo.lt.0.0d0) get_charge_2D(2)=qtopo
get_charge_2D(3:5)=vort

end function get_charge_2D


function get_charge_lat(lat,neigh)result(res)
    use m_vector, only : cross,norm,area
    Implicit none
    real(8)                    :: res(5)
    integer, intent(in)        :: neigh(:,:)
    type(lattice), intent(in)  :: lat
    ! internal
    real(8) :: vort(3),qtopo
    real(8),pointer :: vec(:,:)
    integer     ::  i
    
    res=0.0d0
    vec(1:3,1:product(lat%dim_lat)*lat%nmag)=>lat%M%all_modes

    do i=1,size(neigh,2)
        qtopo=area(vec(:,neigh(1,i)),vec(:,neigh(2,i)),vec(:,neigh(3,i)))+area(vec(:,neigh(1,i)),vec(:,neigh(3,i)),vec(:,neigh(4,i)))
        
        vort=cross(vec(:,neigh(1,i)),vec(:,neigh(2,i)),1,3)+cross(vec(:,neigh(2,i)),vec(:,neigh(3,i)),1,3) &
                +cross(vec(:,neigh(3,i)),vec(:,neigh(4,i)),1,3)+cross(vec(:,neigh(4,i)),vec(:,neigh(1,i)),1,3)
        
        if (qtopo.gt.0.0d0)then
            res(1)=res(1)+qtopo
        else
            res(2)=res(2)+qtopo
        endif
        res(3:5)=res(3:5)+vort
    end do
    nullify(vec)

end function 


function get_charge_lat_i(i,lat,neigh)result(res)
    use m_vector, only : cross,norm,area
    real(8)                    :: res(5)
    integer, intent(in)        :: i
    integer, intent(in)        :: neigh(:,:)
    type(lattice), intent(in)  :: lat
    ! internal
    real(8) :: vort(3),qtopo
    real(8),pointer :: vec(:,:)
    
    res=0.0d0
    vec(1:3,1:product(lat%dim_lat)*lat%nmag)=>lat%M%all_modes
    
    qtopo=area(vec(:,neigh(1,i)),vec(:,neigh(2,i)),vec(:,neigh(3,i)))+area(vec(:,neigh(1,i)),vec(:,neigh(3,i)),vec(:,neigh(4,i)))
    
    vort=cross(vec(:,neigh(1,i)),vec(:,neigh(2,i)),1,3)+cross(vec(:,neigh(2,i)),vec(:,neigh(3,i)),1,3) &
            +cross(vec(:,neigh(3,i)),vec(:,neigh(4,i)),1,3)+cross(vec(:,neigh(4,i)),vec(:,neigh(1,i)),1,3)
    
    if (qtopo.gt.0.0d0)then
        res(1)=qtopo
    else
        res(2)=qtopo
    endif
    res(3:5)=vort
    nullify(vec)
    
end function 


!function get_charge_2D(i,lat)
!use m_vector, only : cross,norm,area
!use m_derived_types, only: lattice
!Implicit none
!real(kind=8), dimension(5) :: get_charge_2D
!integer, intent(in) :: i
!type(lattice),intent(in)    ::  lat
!! internal
!real(kind=8) :: vort(3),qtopo
!
!qtopo=0.0d0
!vort=0.0d0
!get_charge_2D=0.0d0
!
!qtopo=area(Q_topo(1,i)%w,Q_topo(2,i)%w,Q_topo(3,i)%w)+area(Q_topo(1,i)%w,Q_topo(3,i)%w,Q_topo(4,i)%w)
!
!vort=cross(Q_topo(1,i)%w,Q_topo(2,i)%w,1,3)+cross(Q_topo(2,i)%w,Q_topo(3,i)%w,1,3) &
!        +cross(Q_topo(3,i)%w,Q_topo(4,i)%w,1,3)+cross(Q_topo(4,i)%w,Q_topo(1,i)%w,1,3)
!
!if (qtopo.gt.0.0d0) get_charge_2D(1)=qtopo
!if (qtopo.lt.0.0d0) get_charge_2D(2)=qtopo
!get_charge_2D(3:5)=vort
!
!STOP 'get_all_2D'
!
!end function get_charge_2D

!
! subroutine that associates the the topological charge operator and the different vectors
!
subroutine associate_Q_operator(spins,Periodic_log,size_spins)
use m_vector, only : norm
use m_get_position
implicit none
type(vec_point), intent(in) :: spins(:)
logical, intent(in) :: Periodic_log(:)
integer, intent(in) :: size_spins(:)
! internal variables
integer :: i_x,i_y,i_z,i_m
integer :: Ilat(4),position
integer :: Ilat_vois(4),position_vois
integer :: ipu,ipv

do i_m=1,size_spins(4)
Ilat(4)=i_m
Ilat_vois(4)=i_m

   do i_z=1,size_spins(3)
   Ilat(3)=i_z
   Ilat_vois(3)=i_z

      do i_y=1,size_spins(2)
      Ilat(2)=i_y
      select case (size_spins(2))
       case (1)
        ipv=i_y
       case default
        if (Periodic_log(2)) then
         ipv=mod(i_y+1+size_spins(2)-1,size_spins(2))+1
        else
         ipv=i_y
         if (i_y.lt.size_spins(2)) ipv=i_y+1
        endif
      end select

         do i_x=1,size_spins(1)
         Ilat(1)=i_x
         select case (size_spins(1))
          case (1)
           ipu=i_x
          case default
           if (Periodic_log(1)) then
            ipu=mod(i_x-1+1+size_spins(1),size_spins(1))+1
           else
            ipu=i_x
            if (i_x.lt.size_spins(1)) ipu=i_x+1
           endif
         end select

         position=get_position_ND_to_1D(Ilat,size_spins)
         Q_topo(1,position)%w=>spins(position)%w

         ! ipu
         Ilat_vois(1)=ipu
         Ilat_vois(2)=i_y
         position_vois=get_position_ND_to_1D(Ilat_vois,size_spins)
         Q_topo(2,position)%w=>spins(position_vois)%w

         !ipuv
         Ilat_vois(1)=ipu
         Ilat_vois(2)=ipv
         position_vois=get_position_ND_to_1D(Ilat_vois,size_spins)
         Q_topo(3,position)%w=>spins(position_vois)%w

         !ipv
         Ilat_vois(1)=i_x
         Ilat_vois(2)=ipv
         position_vois=get_position_ND_to_1D(Ilat_vois,size_spins)
         Q_topo(4,position)%w=>spins(position_vois)%w


         enddo
      enddo
   enddo
enddo

end subroutine associate_Q_operator

subroutine get_next_ind(lat_size,periodic,ind_next)
    !gets the next index considering periodicity and size=1-> 1
    integer,allocatable ,intent(out)    :: ind_next(:)
    integer,intent(in)                  :: lat_size
    logical,intent(in)                  :: periodic

    integer     :: i

    allocate(ind_next(lat_size),source=1)
    if(lat_size==1) return
    do i=1,lat_size
        ind_next(i)=i+1
    enddo
    if(periodic)then
        ind_next(lat_size)=1
    else
        ind_next(lat_size)=lat_size
    endif
end subroutine


subroutine neighbor_Q(lat,neigh)
    !fill neigh-array with neighboring indices in (3,:) format for magnetization
    !with the neighbors (0,0), (1,0), (1,1), (0,1) for calculating the topological charge
    type(lattice),intent(in)        :: lat
    integer,allocatable,intent(out) :: neigh(:,:)
    ! internal variables
    integer :: i_x,i_y,i_z,i_m
    integer :: pos, Ilat(4)            !position in mag(3,:) and (im,ix,iy,iz)-format
    integer :: pos_vois,Ilat_vois(4)   !position in mag(3,:) and (im,ix,iy,iz)-format
    integer :: size_spins(4)

    integer,allocatable :: vois_x(:),vois_y(:) !neighbors ix/iy index in x and y direction
    
    size_spins(1:3)=lat%dim_lat
    size_spins(4)=lat%M%dim_mode/3 !nmag
    allocate(neigh(4,product(size_spins)),source=0)

    Call get_next_ind(lat%dim_lat(1),lat%periodic(1),vois_x)
    Call get_next_ind(lat%dim_lat(2),lat%periodic(2),vois_y)

    do i_m=1,size_spins(4)
      Ilat(1)=i_m
      Ilat_vois(1)=i_m
      do i_z=1,size_spins(3)
        Ilat(4)=i_z
        Ilat_vois(4)=i_z
        do i_y=1,size_spins(2)
          Ilat(3)=i_y
          do i_x=1,size_spins(1)
            Ilat(2)=i_x
   
            pos=lat%index_4_1(Ilat)
            neigh(1,pos)=pos
    
            ! neighbor (1,0)
            Ilat_vois(2)=vois_x(i_x)
            Ilat_vois(3)=i_y
            pos_vois=lat%index_4_1(Ilat_vois)
            neigh(2,pos)=pos_vois
    
            ! neighbor (1,1)
            Ilat_vois(2)=vois_x(i_x)
            Ilat_vois(3)=vois_y(i_y)
            pos_vois=lat%index_4_1(Ilat_vois)
            neigh(3,pos)=pos_vois
    
            ! neighbor (0,1)
            Ilat_vois(2)=i_x
            Ilat_vois(3)=vois_y(i_y)
            pos_vois=lat%index_4_1(Ilat_vois)
            neigh(4,pos)=pos_vois
          enddo
        enddo
      enddo
    enddo
    deallocate(vois_x,vois_y)
    
end subroutine

!
! get the size of the topological charge operator depending on the geometry
!
subroutine get_size_Q_operator(my_lattice)
implicit none
type(lattice), intent(in) :: my_lattice
!internal
integer :: n_corner,lattice_size(4),N_site
integer :: i,j

lattice_size=shape(my_lattice%ordpar%l_modes)
N_site=product(lattice_size)
n_corner=0

if (.not.my_lattice%periodic(3)) then
   write(6,'(a)') 'no z-periodic boundary conditions'
   n_corner=4*lattice_size(4)

   allocate(Q_topo(n_corner,N_site))

endif

do i=1,N_site
   do j=1,n_corner
      nullify(Q_topo(j,i)%w)
   enddo
enddo

end subroutine get_size_Q_operator

end module m_topo_commons
