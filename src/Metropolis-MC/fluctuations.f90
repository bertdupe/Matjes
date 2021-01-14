module m_fluct
use m_derived_types,only: lattice
use m_get_position
implicit none
private
public fluct_parameters,init_fluct_parameter,eval_fluct

type fluct_parameters
    logical                 :: l_use=.False. 
    integer, allocatable    :: flat_nei(:,:)
    integer, allocatable    :: indexNN(:)
end type
contains

subroutine init_fluct_parameter(fluct_val,lat,l_use_in)
	use m_convert, only : convert
	use m_io_files_utils, only: open_file_write,close_file
    type(fluct_parameters)     :: fluct_val
    type(lattice),intent(in)    :: lat
    logical,intent(in)          :: l_use_in

    fluct_val%l_use=l_use_in
    if(fluct_val%l_use)then
        Call get_neighbours(lat,fluct_val%flat_nei,fluct_val%indexNN)

    else
        !allocate to get more meaningfull error messages in case of mistakes
        allocate(fluct_val%flat_nei(1,1),fluct_val%indexNN(1),source=654123)
    endif

end subroutine



subroutine eval_fluct(MjpMim_sum, MjpMim_ij_sum,  MipMip_sum, MipMim_sum, MipMjp_sum,lat,fluct_val)

    use m_vector, only : cross
    complex(8),intent(inout)           :: MjpMim_sum(:), MjpMim_ij_sum(:,:), MipMip_sum, MipMim_sum,MipMjp_sum
    type(lattice),intent(in)           :: lat
    type(fluct_parameters),intent(in)  :: fluct_val

    integer         :: N_cell
    integer         :: i_cell, i_sh,i_nei,j_flat 
    real(8),pointer :: M3(:,:)
    real(8)         :: Mi(3), Mj(3)
    complex(8),allocatable	 	:: MjpMim(:),  MjpMim_ij(:,:) 
	complex(8)		:: MipMip, MipMim, MipMjp


    if(.not.fluct_val%l_use) return
	i_sh=1
    N_cell=lat%Ncell
    M3(1:3,1:lat%nmag*product(lat%dim_lat))=>lat%M%all_modes
    if(lat%nmag>1) STOP "FLUCTATIONS WILL NOT WORK WITH MORE THAN ONE MAGNETIC MOMENT IN THE UNIT_CELL"
	allocate(MjpMim(fluct_val%indexNN(i_sh)))
	allocate(MjpMim_ij(fluct_val%indexNN(i_sh),N_cell)) !store average per unique pair of neighbours

    MjpMim(:)=cmplx(0.0d0,0.0d0,8);    MipMip=cmplx(0.0d0,0.0d0,8)
    MipMim=cmplx(0.0d0,0.0d0,8);    MipMjp=cmplx(0.0d0,0.0d0,8)
	MjpMim_ij(:,:)=cmplx(0.0d0,0.0d0,8)

    do i_cell=1,N_cell  !loop on sites
        Mi=M3(:,i_cell) !site i

        MipMip=MipMip+cmplx(Mi(1)**2-Mi(2)**2, 2.0d0*Mi(1)*Mi(2),8)
        MipMim=MipMim+cmplx(Mi(1)**2+Mi(2)**2,     0.0d0        ,8)

        do i_sh=1,size(fluct_val%indexNN) !loop on shell of neighbours
            do i_nei=1,fluct_val%indexNN(i_sh) !loop on neighbours
                if(fluct_val%flat_nei(i_cell,i_nei).eq.-1) cycle !skip non-connected neighbours
                j_flat=fluct_val%flat_nei(i_cell,i_nei)
                Mj=M3(:,j_flat) !site j


				MjpMim_ij(i_nei,i_cell)=cmplx(Mi(1)*Mj(1)+Mi(2)*Mj(2)  ,  Mi(1)*Mj(2)-Mi(2)*Mj(1),8)
                MjpMim(i_nei)=MjpMim(i_nei)+MjpMim_ij(i_nei,i_cell)  ! Mj+Mi-


                MipMjp=MipMjp+cmplx(Mi(1)*Mj(1)-Mi(2)*Mj(2)  ,  Mi(1)*Mj(2)+Mi(2)*Mj(1),8)  ! Mi+Mj+
            enddo
        enddo
    enddo
    MjpMim_sum(:)=MjpMim_sum(:)+MjpMim(:)
    MjpMim_ij_sum(:,:)=MjpMim_ij_sum(:,:)+MjpMim_ij(:,:)
    MipMip_sum=MipMip_sum+MipMip
    MipMim_sum=MipMim_sum+MipMim
    MipMjp_sum=MipMjp_sum+MipMjp

end subroutine


!to allocate fluctuation arrays 
!subroutine allocate_fluct(this,lat,fluct_val)
!	type(lattice),intent(in)           :: lat
 !   type(fluct_parameters),intent(in)  :: fluct_val

!	allocate(this%MjpMim_ij_sum(fluct_val%indexNN(i_sh),lat%Ncell)) 
!	allocate(this%MjpMim_ij_av(fluct_val%indexNN(i_sh),lat%Ncell)) !fluct average per unique pair of neighbours, shape is Nnei x N

!end subroutine allocate_fluct



subroutine get_neighbours(lat,flat_nei, indexNN)
	use m_io_utils
	use m_io_files_utils

    use m_get_table_nn,only :get_table_nn
    Implicit none
    type(lattice),intent(in)  :: lat
    integer,allocatable,intent(inout)        :: flat_nei(:,:), indexNN(:)
    ! internal variables
    integer      :: N_shells,N_cell
    integer, allocatable :: tableNN(:,:,:,:,:,:)
    integer      :: i_sh, i_x, i_y, i_z, i_nei, temp1(3),temp2(3), i_flat, j_flat
    real(kind=8), allocatable :: position(:,:)
	integer :: N
	integer :: propag_dir=1 !direction of propagation of fluctuations to look at, 1,2,3 for x,y,z

    N_cell=lat%Ncell
    N_shells = 1 !choose how many shells of neighbours: 1 = first neighbours  /!\currently boundary wont work for more than 1 shell 
    Call get_table_nn(lat,N_shells,indexNN,tableNN)
    allocate(flat_nei(N_cell,sum(indexNN))) !flat nei is N x number of neighbours per site
	flat_nei(:,:)=0;
	

	!get positions for no reason other than to know in which spatial direction i move for pair of neighbours
	N=get_lines('positions.dat')
	allocate(position(3,N))
	position=0.0d0
	call get_position(position,'positions.dat')
	

	!write(*,*) 'Only fluctuations propagating in lattice direction ',propag_dir,' will be considered.'
    !convert to a table with linear indices (Ncell, indexNN), /!\ taking only forward neighbours along x,y,z
	!simple cubic lattice: columns 1:6 of flat_nei correspond to neighbours along +-x,+-y,+-z
    do i_sh=1,N_shells !loop on shell of neighbours
        do i_x=1,lat%dim_lat(1) !loop on lattice sites
            do i_y=1,lat%dim_lat(2)
                do i_z=1,lat%dim_lat(3) 
                    temp1=[i_x,i_y,i_z]
                    i_flat=lat%index_m_1(temp1)  !get i_flat
					!write(*,*) 'spin: ',i_flat,' with xyz positions ', position(:,i_flat) 
                    do i_nei=1,indexNN(i_sh) !loop on neighbours
						
						!write(*,*)'spin',i_flat, ' with  coordinates ',temp1(1),temp1(2),temp1(3) 
                        temp2=tableNN(1:3,i_nei,i_x,i_y,i_z,1)!get x,y,z indices of neighbour


						if( (all(temp1.ge.temp2)).or.(tableNN(5,i_nei,i_x,i_y,i_z,1).ne.1) ) then  !do not keep backwards neighbours or non-connected neighbours
							j_flat=-1; 
						else if ( ( temp1(1)==1.and.temp2(1)==lat%dim_lat(1).and.lat%dim_lat(1).ne.1 ) &
						.or.      ( temp1(2)==1.and.temp2(2)==lat%dim_lat(2).and.lat%dim_lat(2).ne.1 ) &
						.or.      ( temp1(3)==1.and.temp2(3)==lat%dim_lat(3).and.lat%dim_lat(3).ne.1 ) )  then !boundary:remove left and top
							j_flat=-1 
							!write(*,*) 'boundary nei to remove: verified for temp1, temp2=', temp1(1:3), temp2(1:3)
						else
                       		j_flat=lat%index_m_1(temp2) !get j_flat
						end if

						!boundary: add right and bottom
						if ( (temp1(1)==lat%dim_lat(1).and.temp2(1)==1.and.lat%dim_lat(1).ne.1) &
						.or. (temp1(2)==lat%dim_lat(2).and.temp2(2)==1.and.lat%dim_lat(2).ne.1) &
						.or. (temp1(3)==lat%dim_lat(3).and.temp2(3)==1.and.lat%dim_lat(3).ne.1)	) then
                       		j_flat=lat%index_m_1(temp2) 
						end if

						!only keep neighbour along the propagation direction
						!if (temp1(propag_dir).eq.temp2(propag_dir)) then
						!	j_flat=-1 
						!endif
					

                        flat_nei(i_flat,i_nei)=j_flat

						!write(*,*)'>>neighbour',j_flat, ' with xyz positions ', position(:,j_flat) 
                    enddo
                enddo
            enddo
        enddo
    enddo

!do i_flat=1,N_cell
!	write(*,*)'flat_nei for spin',i_flat, ' is', flat_nei(i_flat,:)
!enddo

END subroutine get_neighbours
end module 

