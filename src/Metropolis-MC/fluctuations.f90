module m_fluct
use m_derived_types,only: lattice
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



subroutine eval_fluct(MjpMim_sum, MipMip_sum, MipMim_sum, MipMjp_sum,lat,fluct_val)

    use m_vector, only : cross
    complex(8),intent(inout)           :: MjpMim_sum, MipMip_sum, MipMim_sum,MipMjp_sum
    type(lattice),intent(in)           :: lat
    type(fluct_parameters),intent(in)  :: fluct_val

    real(8),pointer :: M3(:,:)
    real(8)         :: Mi(3), Mj(3)
    complex(8)      :: MjpMim, MipMip, MipMim, MipMjp
    integer         :: N_cell
    integer         :: i_cell, i_sh,i_nei,j_flat 

    if(.not.fluct_val%l_use) return
    N_cell=lat%Ncell
    M3(1:3,1:lat%nmag*product(lat%dim_lat))=>lat%M%all_modes
    if(lat%nmag>1) STOP "FLUCTATIONS WILL NOT WORK WITH MORE THAN ONE MAGNETIC MOMENT IN THE UNIT_CELL"

    MjpMim=cmplx(0.0d0,0.0d0,8);    MipMip=cmplx(0.0d0,0.0d0,8)
    MipMim=cmplx(0.0d0,0.0d0,8);    MipMjp=cmplx(0.0d0,0.0d0,8)

    do i_cell=1,N_cell  !loop on sites
        Mi=M3(:,i_cell) !site i

        MipMip=MipMip+cmplx(Mi(1)**2-Mi(2)**2, 2.0d0*Mi(1)*Mi(2),8)
        MipMim=MipMim+cmplx(Mi(1)**2+Mi(2)**2,     0.0d0        ,8)

        do i_sh=1,size(fluct_val%indexNN) !loop on shell of neighbours
            do i_nei=1,fluct_val%indexNN(i_sh) !loop on neighbours

                if(fluct_val%flat_nei(i_cell,i_nei).eq.-1) cycle !skip non-connected neighbours
                j_flat=fluct_val%flat_nei(i_cell,i_nei)
                Mj=M3(:,j_flat) !site j

                MjpMim=MjpMim+cmplx(Mi(1)*Mj(1)+Mi(2)*Mj(2)  ,  Mi(1)*Mj(2)-Mi(2)*Mj(1),8)  ! Mj+Mi-
                MipMjp=MipMjp+cmplx(Mi(1)*Mj(1)-Mi(2)*Mj(2)  ,  Mi(1)*Mj(2)+Mi(2)*Mj(1),8)  ! Mi+Mj+


				!write(*,*) 'i,j  = ', i_cell, j_flat,  ', Mi = ' , Mi(1), Mi(2),Mi(3), 'Mj=', Mj(1), Mj(2),Mj(3), '  and (Mi x Mj)_z =  ', aimag( MjpMim)
            enddo
        enddo

    enddo
    MjpMim_sum=MjpMim_sum+MjpMim
    MipMip_sum=MipMip_sum+MipMip
    MipMim_sum=MipMim_sum+MipMim
    MipMjp_sum=MipMjp_sum+MipMjp

	!write(*,*) 'end of one step, MjpMim_sum = ' , MjpMim_sum

end subroutine

subroutine get_neighbours(lat,flat_nei, indexNN)
    use m_get_table_nn,only :get_table_nn
    Implicit none
    type(lattice),intent(in)  :: lat
    integer,allocatable,intent(inout)        :: flat_nei(:,:), indexNN(:)
    ! internal variables
    integer      :: N_shells,N_cell
    integer, allocatable :: tableNN(:,:,:,:,:,:)
    integer      :: i_sh, i_x, i_y, i_z, i_nei, temp1(3),temp2(3), i_flat, j_flat
    
    N_cell=lat%Ncell
    N_shells = 1 !choose how many shells of neighbours: 1 = first neighbours
    Call get_table_nn(lat,N_shells,indexNN,tableNN)
    allocate(flat_nei(N_cell,sum(indexNN))) !flat nei is N x number of neighbours per site
	flat_nei(:,:)=0;

    !convert to a table with linear indices (Ncell, indexNN), /!\ taking only forward neighbours along x,y,z
    do i_sh=1,N_shells !loop on shell of neighbours
        do i_x=1,lat%dim_lat(1) !loop on lattice sites
            do i_y=1,lat%dim_lat(2)
                do i_z=1,lat%dim_lat(3) 
                    temp1=[i_x,i_y,i_z]
                    i_flat=lat%index_m_1(temp1)  !get i_flat
                    do i_nei=1,indexNN(i_sh) !loop on neighbours
                        temp2=tableNN(1:3,i_nei,i_x,i_y,i_z,1)!get x,y,z indices of neighbour
				
						!write(*,*) 'temp1=', temp1(1:2), ' temp2= ', temp2(1:2), ' temp1.gt.temp2 = ' , all(temp1.ge.temp2)

						if( (all(temp1.ge.temp2)).or.(tableNN(5,i_nei,i_x,i_y,i_z,1).ne.1) ) then  !do not keep backwards neighbours or non-connected neighbours
							j_flat=-1; 
						else if ( ( temp1(2)==1.and.temp2(2)==lat%dim_lat(2) ) .or. ( temp1(1)==1.and.temp2(1)==lat%dim_lat(1)) ) then !boundary:remove left and top
							j_flat=-1 
						else
                       		j_flat=lat%index_m_1(temp2) !get j_flat
						end if

						!boundary: add right and bottom
						if ( (temp1(2)==lat%dim_lat(2).and.temp2(2)==1) .or. (temp1(1)==lat%dim_lat(1).and.temp2(1)==1) ) then
                       		j_flat=lat%index_m_1(temp2) 
						end if

                        flat_nei(i_flat,i_nei)=j_flat
                    enddo
                enddo
            enddo
        enddo
    enddo

! do i_flat=1,N_cell
!	write(*,*)'flat_nei for spin',i_flat, ' =', flat_nei(i_flat,:)
!enddo

END subroutine get_neighbours
end module 

