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


! ------- for <M+M->
subroutine eval_fluct(Re_MpMm_sum,Im_MpMm_sum,lat,fluct_val)
    use m_vector, only : cross
    real(8),intent(inout)               :: Re_MpMm_sum,Im_MpMm_sum
    type(lattice),intent(in)            :: lat
    type(fluct_parameters),intent(in)  :: fluct_val

    real(8),pointer :: M3(:,:)
    real(8)         :: temp(3)
    real(8)         :: Mi(3), Mj(3), Re_MpMm, Im_MpMm
    integer         :: N_neighbours,N_cell
    integer         :: i_cell, i_sh,i_nei,j_flat 

    if(.not.fluct_val%l_use) return
	N_neighbours=rank(fluct_val%indexNN) !number of shells of neighbours to be used
    N_cell=lat%Ncell
    M3(1:3,1:lat%nmag*product(lat%dim_lat))=>lat%M%all_modes
    Re_MpMm=0.0d0
    Im_MpMm=0.0d0
    do i_cell=1,N_cell  !loop on sites
        Mi=M3(:,i_cell) !site i
        do i_sh=1,N_neighbours !loop on shell of neighbours
            do i_nei=1,fluct_val%indexNN(i_sh) !loop on neighbours
                if(fluct_val%flat_nei(i_cell,i_nei).eq.-1) cycle !skip non-connected neighbours
                j_flat=fluct_val%flat_nei(i_cell,i_nei)

                !write(*,*) "site", i_cell, "neighbour", j_flat
                Mj=M3(:,j_flat) !site j
                Re_MpMm = Re_MpMm + dot_product(Mi,Mj) - Mi(3)*Mj(3) !real part M+M-, sum over all pairs
                temp = cross(Mi,Mj)
                Im_MpMm= Im_MpMm + temp(3) !imaginary part of M+M-, sum over all pairs
				!write(*,*)'Mi=',Mi(1:3), ' Mj= ',Mj(1:3), 'Mi*Mj=',dot_product(Mi,Mj) ,' Mi*Mj-MizMjz = ',dot_product(Mi,Mj) - Mi(3)*Mj(3),' (MixMj)_z = ', temp(3)
				!write(*,*) 'Im_MpMm = ', Im_MpMm
            enddo
        enddo
    enddo
    Re_MpMm_sum=Re_MpMm_sum+Re_MpMm
    Im_MpMm_sum=Im_MpMm_sum+Im_MpMm
	!write(*,*) 'in eval_fluct, Re_MpMm_sum= ',Re_MpMm_sum,' Im_MpMm_sum= ',Im_MpMm_sum

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
    !write(*,*) "N_cell= " ,N_cell, "sum(indexNN)= ", sum(indexNN)
    allocate(flat_nei(N_cell,N_cell*sum(indexNN))) !flat nei is N x N*number of neighbours per site
    !write(*,*) "in get_nei l223 shape(flat_nei)=", shape(flat_nei) 

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
                        !write(*,*) i_sh, i_x, i_y, i_z, temp(1:3), i_flat, i_nei, j_flat
                        !write(*,*) "shape(flat_nei)=", shape(flat_nei)
                        flat_nei(i_flat,i_nei)=j_flat
                    enddo
                enddo
            enddo
        enddo
    enddo
        

END subroutine get_neighbours
end module 

