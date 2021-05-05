module m_fluct
use m_derived_types,only: lattice
!use m_get_position
use m_neighbor_type
implicit none
private
public fluct_parameters,init_fluct_parameter,eval_fluct

type fluct_parameters
    logical                 :: l_use=.False. 
!    integer, allocatable    :: flat_nei(:,:)
!    integer, allocatable    :: indexNN(:)
    type(neighbors)         :: neigh
    integer                 :: ind(1)  !index of positive and negative nearest neighbor in neigh%diff_vec array 
contains
    procedure get_nneigh
end type
contains

function get_nneigh(this)result(nneigh)
    class(fluct_parameters),intent(in)  :: this
    integer                             :: nneigh

    nneigh=this%neigh%Nshell(2) !2 because one is on-site, 2 next nearest neighbor as decided in init_fluct_parameter
end function

subroutine init_fluct_parameter(fluct_val,lat,l_use_in,direction)
    use, intrinsic :: iso_fortran_env, only : output_unit
    type(fluct_parameters)     :: fluct_val
    type(lattice),intent(in)    :: lat
    logical,intent(in)          :: l_use_in
    real(8),intent(in)          :: direction(3)

    integer                     :: bnd(2)
    real(8),allocatable         :: diff_vec(:,:)
    real(8),allocatable         :: proj(:) 

    fluct_val%l_use=l_use_in
    if(fluct_val%l_use)then
        Call fluct_val%neigh%get([1,1],[0,1],lat)   !so far only fluctuations implemented between next nearest neighbors of first atom type
        bnd=[fluct_val%neigh%Nshell(1)+1,sum(fluct_val%neigh%Nshell(1:2))]
        allocate(diff_vec(3,bnd(2)-bnd(1)+1))
        diff_vec=fluct_val%neigh%diff_vec(:,bnd(1):bnd(2))
        allocate(proj(size(diff_vec,2)),source=0.0d0)
        proj=matmul(direction,diff_vec)
        fluct_val%ind(1)=maxloc(proj,1)
       ! proj=matmul(diff_vec(:,fluct_val%ind(1)),diff_vec)
        !fluct_val%ind(2)=minloc(proj,1)
        fluct_val%ind=fluct_val%ind+bnd(1)-1

        write(output_unit,'(//A)') "Fluctuation neighbors are in the following order:"
        write(output_unit,'(3F16.8)') diff_vec
        write(output_unit,'(/A)') "Considering the following neighbors:"
        write(output_unit,'(3F16.8)') fluct_val%neigh%diff_vec(:,fluct_val%ind)

    endif
end subroutine

subroutine eval_fluct( MjpMim_ij_sum, MipMip_sum, MipMim_sum, MipMjp_sum,lat,fluct_val)
    use m_vector, only : cross
    complex(8),intent(inout)           :: MjpMim_ij_sum(:,:), MipMip_sum, MipMim_sum,MipMjp_sum
    type(lattice),intent(in)           :: lat
    type(fluct_parameters),intent(in)  :: fluct_val

    integer         :: N_cell
    integer         :: i_cell, i_sh,i_nei,j_flat 
    integer         :: i,j
    real(8)         :: Mi(3), Mj(3)
	complex(8)		:: MipMip, MipMim, MipMjp
    integer         :: neigh_start
    integer         :: shell_offset
    integer         :: ind

    if(.not.fluct_val%l_use) return

    MipMip=cmplx(0.0d0,0.0d0,8)
    MipMjp=cmplx(0.0d0,0.0d0,8)
    MipMim=cmplx(0.0d0,0.0d0,8)

    !do onsite-terms
    neigh_start=1
    do i_sh=1,fluct_val%neigh%Nshell(1)
        do i_nei=neigh_start,fluct_val%neigh%ishell(i_sh) 
            i=fluct_val%neigh%pairs(1,i_nei)
            Mi=lat%M%modes_3(:,i)
            MipMip=MipMip+cmplx(Mi(1)**2-Mi(2)**2, 2.0d0*Mi(1)*Mi(2),8)
            MipMim=MipMim+cmplx(Mi(1)**2+Mi(2)**2,     0.0d0        ,8)
        enddo
        neigh_start=fluct_val%neigh%ishell(i_sh)+1
    enddo

    !do nearest neighbor terms
    do i_sh=1,size(fluct_val%ind)
        ind=fluct_val%ind(i_sh)
        do i_nei=fluct_val%neigh%ishell(ind-1)+1, fluct_val%neigh%ishell(ind) 
            i=fluct_val%neigh%pairs(1,i_nei)
            j=fluct_val%neigh%pairs(2,i_nei)
            Mi=lat%M%modes_3(:,i)
            Mj=lat%M%modes_3(:,j)
            MipMjp               =MipMjp               +cmplx(Mi(1)*Mj(1)-Mi(2)*Mj(2)  ,  Mi(1)*Mj(2)+Mi(2)*Mj(1),8)     !why is here a minus?
            MjpMim_ij_sum(i_sh,i)=MjpMim_ij_sum(i_sh,i)+cmplx(Mi(1)*Mj(1)+Mi(2)*Mj(2)  ,  Mi(1)*Mj(2)-Mi(2)*Mj(1),8)
        enddo
    enddo

    MipMip_sum=MipMip_sum+MipMip
    MipMim_sum=MipMim_sum+MipMim
    MipMjp_sum=MipMjp_sum+MipMjp
end subroutine


!subroutine get_neighbours(lat,flat_nei, indexNN)
!	use m_io_utils
!	use m_io_files_utils
!
!    use m_get_table_nn,only :get_table_nn
!    Implicit none
!    type(lattice),intent(in)  :: lat
!    integer,allocatable,intent(inout)        :: flat_nei(:,:), indexNN(:)
!    ! internal variables
!    integer      :: N_shells,N_cell
!    integer, allocatable :: tableNN(:,:,:,:,:,:)
!    integer      :: i_sh, i_x, i_y, i_z, i_nei, temp1(3),temp2(3), i_flat, j_flat
!    real(kind=8), allocatable :: position(:,:)
!	integer :: N
!	integer :: propag_dir=1 !direction of propagation of fluctuations to look at, 1,2,3 for x,y,z
!
!    N_cell=lat%Ncell
!    N_shells = 1 !choose how many shells of neighbours: 1 = first neighbours  /!\currently boundary wont work for more than 1 shell 
!    Call get_table_nn(lat,N_shells,indexNN,tableNN)
!    allocate(flat_nei(N_cell,sum(indexNN))) !flat nei is N x number of neighbours per site
!    flat_nei=0
!
!
!	!get positions for no reason other than to know in which spatial direction the bond is for given pair of neighbours
!	N=get_lines('positions.dat')
!	allocate(position(3,N))
!	position=0.0d0
!	call get_position(position,'positions.dat')
!	
!
!	!simple cubic lattice: columns 1:6 of flat_nei(i,:) correspond to neighbours of i along: x,y,z,-z,-y,-x
!	!careful: if film in xy plane then it's only along x,y,-y,-x
!    do i_sh=1,N_shells !loop on shell of neighbours
!        do i_x=1,lat%dim_lat(1) !loop on lattice sites
!            do i_y=1,lat%dim_lat(2)
!                do i_z=1,lat%dim_lat(3) 
!                    temp1=[i_x,i_y,i_z]
!                    i_flat=lat%index_m_1(temp1)  !get i_flat
!					!write(*,*) 'spin: ',i_flat,' with xyz positions ', position(:,i_flat) 
!                    do i_nei=1,indexNN(i_sh) !loop on neighbours
!						
!						!write(*,*)'spin',i_flat, ' with  coordinates ',temp1(1),temp1(2),temp1(3) 
!                        temp2=tableNN(1:3,i_nei,i_x,i_y,i_z,1)!get x,y,z indices of neighbour
!
!						if (tableNN(5,i_nei,i_x,i_y,i_z,1).ne.1) j_flat=-1 !remove non-connected neighbours
!						!if( (all(temp1.ge.temp2)).or.(tableNN(5,i_nei,i_x,i_y,i_z,1).ne.1) ) then  !do not keep backwards neighbours or non-connected neighbours
!						!	j_flat=-1; 
!						!else if ( ( temp1(1)==1.and.temp2(1)==lat%dim_lat(1).and.lat%dim_lat(1).ne.1 ) &
!						!.or.      ( temp1(2)==1.and.temp2(2)==lat%dim_lat(2).and.lat%dim_lat(2).ne.1 ) &
!						!.or.      ( temp1(3)==1.and.temp2(3)==lat%dim_lat(3).and.lat%dim_lat(3).ne.1 ) )  then !boundary:remove left and top
!						!	j_flat=-1 
!							!write(*,*) 'boundary nei to remove: verified for temp1, temp2=', temp1(1:3), temp2(1:3)
!						!else
!                       		j_flat=lat%index_m_1(temp2) !get j_flat
!						!end if
!
!						!boundary: add right and bottom
!						!if ( (temp1(1)==lat%dim_lat(1).and.temp2(1)==1.and.lat%dim_lat(1).ne.1) &
!						!.or. (temp1(2)==lat%dim_lat(2).and.temp2(2)==1.and.lat%dim_lat(2).ne.1) &
!						!.or. (temp1(3)==lat%dim_lat(3).and.temp2(3)==1.and.lat%dim_lat(3).ne.1)	) then
!                       !		j_flat=lat%index_m_1(temp2) 
!						!end if
!
!						!only keep neighbour along the propagation direction
!						!if (temp1(propag_dir).eq.temp2(propag_dir)) then
!						!	j_flat=-1 
!						!endif
!					
!						if (i_flat.eq.j_flat) j_flat=-1 !remove self interaction
!                        flat_nei(i_flat,i_nei)=j_flat
!
!					!	write(*,*)'>>neighbour',j_flat, ' with xyz positions ', position(:,j_flat) 
!                    enddo
!                enddo
!            enddo
!        enddo
!    enddo
!
!!do i_flat=1,N_cell
!!	write(*,*)'flat_nei for spin',i_flat, ' is', flat_nei(i_flat,:)
!!enddo
!
!END subroutine get_neighbours
end module 

