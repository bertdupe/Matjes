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



subroutine eval_fluct(Re_MipMjm_sum,Im_MipMjm_sum, &
				    & Re_MipMip_sum,Im_MipMip_sum, &
					& Re_MipMim_sum,			   &
					& Re_MipMjp_sum,Im_MipMjp_sum, &																								
					& lat,fluct_val)

    use m_vector, only : cross
    real(8),intent(inout)              :: Re_MipMjm_sum,Im_MipMjm_sum, Re_MipMip_sum,Im_MipMip_sum,Re_MipMim_sum,Re_MipMjp_sum,Im_MipMjp_sum
    type(lattice),intent(in)           :: lat
    type(fluct_parameters),intent(in)  :: fluct_val

    real(8),pointer :: M3(:,:)
    real(8)         :: temp(3)
    real(8)         :: Mi(3), Mj(3), Re_MipMjm, Im_MipMjm, Re_MipMip, Im_MipMip, Re_MipMim, Re_MipMjp, Im_MipMjp
    integer         :: N_neighbours,N_cell
    integer         :: i_cell, i_sh,i_nei,j_flat 

    if(.not.fluct_val%l_use) return
	N_neighbours=rank(fluct_val%indexNN) !number of shells of neighbours to be used
    N_cell=lat%Ncell
    M3(1:3,1:lat%nmag*product(lat%dim_lat))=>lat%M%all_modes

    Re_MipMjm=0.0d0
    Im_MipMjm=0.0d0
	Re_MipMip=0.0d0
	Im_MipMip=0.0d0
	Re_MipMim=0.0d0	
	Re_MipMjp=0.0d0
	Im_MipMjp=0.0d0

    do i_cell=1,N_cell  !loop on sites
        Mi=M3(:,i_cell) !site i

		Re_MipMip=Re_MipMip + Mi(1)**2-Mi(2)**2
		Im_MipMip_sum=Im_MipMip+ 2*Mi(1)*Mi(2)

		Re_MipMim=Re_MipMim + Mi(1)**2+Mi(2)**2

        do i_sh=1,N_neighbours !loop on shell of neighbours
            do i_nei=1,fluct_val%indexNN(i_sh) !loop on neighbours

                if(fluct_val%flat_nei(i_cell,i_nei).eq.-1) cycle !skip non-connected neighbours
                j_flat=fluct_val%flat_nei(i_cell,i_nei)
                Mj=M3(:,j_flat) !site j

                Re_MipMjm = Re_MipMjm + dot_product(Mi,Mj) - Mi(3)*Mj(3) 
                temp = cross(Mi,Mj)
                Im_MipMjm= Im_MipMjm + temp(3) 

				Re_MipMjp=Re_MipMjp + Mi(1)*Mj(1) -  Mi(2)*Mj(2)
				Im_MipMjp=Im_MipMjp + Mi(1)*Mj(2) +  Mi(2)*Mj(1)

            enddo
        enddo

    enddo
    Re_MipMjm_sum=Re_MipMjm_sum+Re_MipMjm
    Im_MipMjm_sum=Im_MipMjm_sum+Im_MipMjm

	Re_MipMip_sum=Re_MipMip_sum+Re_MipMip
	Im_MipMip_sum=Im_MipMip_sum+Im_MipMip

	Re_MipMim_sum=Re_MipMim_sum+Re_MipMim

	Re_MipMjp_sum=Re_MipMjp_sum+Re_MipMjp
	Im_MipMjp_sum=Im_MipMjp_sum+Im_MipMjp

end subroutine

subroutine get_neighbours(lat,flat_nei, indexNN)
    use m_get_table_nn,only :get_table_nn
    Implicit none
    type(lattice),intent(in)  :: lat
    integer,allocatable,intent(inout)        :: flat_nei(:,:), indexNN(:)
    ! internal variables
    integer      :: N_neighbours,N_cell
    integer, allocatable :: tableNN(:,:,:,:,:,:)
    integer      :: i_sh, i_x, i_y, i_z, i_nei, temp(3), i_flat, j_flat
    
    N_cell=lat%Ncell

    N_neighbours = 1 !choose how many shells of neighbours: 1 = first neighbours
    Call get_table_nn(lat,N_neighbours,indexNN,tableNN)
    allocate(flat_nei(N_cell,N_cell*sum(indexNN))) !flat nei is N x N*number of neighbours per site


    !convert to a table with linear indices (Ncell, indexNN)
    do i_sh=1,N_neighbours !loop on shell of neighbours
        do i_x=1,lat%dim_lat(1) !loop on lattice sites
            do i_y=1,lat%dim_lat(2)
                do i_z=1,lat%dim_lat(3) 
                    temp=[i_x,i_y,i_z]
                    i_flat=lat%index_m_1(temp)  !get i_flat
                    do i_nei=1,indexNN(i_sh) !loop on neighbours
                        temp=tableNN(1:3,i_nei,i_x,i_y,i_z,1)!get x,y,z indices of neighbour
                        j_flat=lat%index_m_1(temp) !get j_flat
                        if(tableNN(5,i_nei,i_x,i_y,i_z,1).ne.1) j_flat=-1 !if neighbours are not connected set to -1
                        flat_nei(i_flat,i_nei)=j_flat
                    enddo
                enddo
            enddo
        enddo
    enddo
        

END subroutine get_neighbours
end module 

