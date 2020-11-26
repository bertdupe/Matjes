module m_topo_commons
use m_derived_types, only : t_cell,lattice
Implicit none

interface get_charge
  module procedure :: get_charge_lat
  module procedure :: get_charge_lat_i
!  module procedure :: get_charge_res_lat

end interface

private
public :: get_charge,neighbor_Q

contains
!!cannot do that since the interface is ambiguous, using subroutines would allow that, but one has to change a few things...
!function get_charge_res_lat(lat,neigh)result(topo_out)
!    use m_vector, only : cross,norm,area
!    Implicit none
!    integer, intent(in)        :: neigh(:,:)
!    type(lattice), intent(in)  :: lat
!    real(8)                    :: topo_out(5,size(neigh,2))
!    ! internal
!    real(8) :: vort(3),qtopo
!    real(8),pointer :: vec(:,:)
!    integer     ::  i
!    
!    vec(1:3,1:product(lat%dim_lat)*lat%nmag)=>lat%M%all_modes
!    do i=1,size(neigh,2)
!        qtopo=area(vec(:,neigh(1,i)),vec(:,neigh(2,i)),vec(:,neigh(3,i)))+area(vec(:,neigh(1,i)),vec(:,neigh(3,i)),vec(:,neigh(4,i)))
!        
!        topo_out(3:5,i)=cross(vec(:,neigh(1,i)),vec(:,neigh(2,i)),1,3)+cross(vec(:,neigh(2,i)),vec(:,neigh(3,i)),1,3) &
!                +cross(vec(:,neigh(3,i)),vec(:,neigh(4,i)),1,3)+cross(vec(:,neigh(4,i)),vec(:,neigh(1,i)),1,3)
!        
!        if (qtopo.gt.0.0d0)then
!            topo_out(1,i)=qtopo
!        else
!            topo_out(2,i)=qtopo
!        endif
!    end do
!    nullify(vec)
!end function

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
    !fill neigh-array with neighboring indices in (4,:) format for magnetization
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

end module m_topo_commons
