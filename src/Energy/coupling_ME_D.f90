module m_coupling_ME_D
use m_input_H_types, only: io_H_ME_D
use m_io_utils, only: get_parameter,get_coeff,number_nonzero_coeff,max_ind_variable
implicit none

private
public :: get_coupling_ME_D, read_ME_D_input

contains

subroutine read_ME_D_input(io_param,fname,io)
    integer,intent(in)              :: io_param
    character(len=*), intent(in)    :: fname
    type(io_H_ME_D),intent(out)     :: io

    Call get_parameter(io_param,fname,'ME_antisym',io%pair,io%is_set) 
end subroutine

subroutine get_coupling_ME_D(Ham,io,lat)
    !get coupling in t_H Hamiltonian format
    use m_H_public
    use m_derived_types
    use m_setH_util,only: get_coo,ind
    use m_neighbor_type, only: neighbors
    use m_mode_public

    class(t_H),intent(inout)    :: Ham
    type(io_H_ME_D),intent(in)  :: io
    type(lattice),intent(in)    :: lat

    !local Hamiltonian
    real(8),allocatable  :: Htmp(:,:)
    !local Hamiltonian in coo format
    real(8),allocatable  :: val_tmp(:)
    integer,allocatable  :: ind_tmp(:,:)

    class(t_H),allocatable    :: Ham_tmp    !temporary Hamiltonian type used to add up Ham

    integer         :: i_atpair,N_atpair    !loop parameters which atom-type connection are considered (different neighbor types)
    integer         :: i_dist,N_dist        !loop parameters which  connection are considered (different neighbor types)
    integer         :: i_pair           !loop keeping track which unique connection between the same atom types is considered (indexes "number shells" in neighbors-type)
    integer         :: i_shell          !counting the number of unique connection for given atom types and a distance
    integer         :: connect_bnd(2)   !indices keeping track of which pairs are used for the particular connection
    type(neighbors) :: neigh            !all neighbor information for a given atom-type pair
    real(8)         :: Hmag             !magnitude of Hamiltonian parameter
    integer         :: atind_mag(2)     !index of considered atom in basis of magnetic atoms (1:Nmag)A
    integer         :: offset_mag(2)    !offset for start in dim_mode of chosed magnetic atom
    integer         :: dim_modes_r(2)   !mode dimension on the right side (M,E)

    !parameters to get the difference between the considered atoms
    integer         :: ind1(4),ind2(4)  !indices to get position distance ([ix,iy,iz,ia]-lattice coordinates)
    real(8)         :: diff_pos(3)  !normalized difference vector between both considered atoms

    if(io%is_set)then
        if(lat%E%dim_mode==0) STOP "E-field has to be set when using coupling_ME"
        Call get_Htype(Ham_tmp)
        N_atpair=size(io%pair)
        dim_modes_r=[lat%M%dim_mode,lat%E%dim_mode]
        allocate(Htmp(lat%M%dim_mode,product(dim_modes_r)))!local Hamiltonian modified for each shell/neighbor
        do i_atpair=1,N_atpair
            !loop over different connected atom types
            Call neigh%get(io%pair(i_atpair)%attype,io%pair(i_atpair)%dist,lat)
            N_dist=size(io%pair(i_atpair)%dist)
            i_pair=0
            connect_bnd=1 !initialization for lower bound
            do i_dist=1,N_dist
                !loop over distances (nearest, next nearest,... neighbor)
                Hmag=io%pair(i_atpair)%val(i_dist)
                do i_shell=1,neigh%Nshell(i_dist)
                    !loop over all different connections with the same distance
                    i_pair=i_pair+1
                    !get the distance vector (diff_pos)
                    ind1(1:3)=lat%index_1_3(neigh%pairs(1,connect_bnd(1)))
                    ind1(4)=neigh%at_pair(1,i_pair)
                    ind2(1:3)=lat%index_1_3(neigh%pairs(2,connect_bnd(1)))
                    ind2(4)=neigh%at_pair(2,i_pair)
                    diff_pos=lat%pos_diff_ind(ind1,ind2) 
                    where(abs(diff_pos)<norm2(diff_pos)*1.0d-8) diff_pos=0.0d0
                    diff_pos=diff_pos/norm2(diff_pos)

                    !set local Hamiltonian
                    !use lagrange identity to write as
                    !(mi.Ei)(mj.r)-(mj.Ei)(mi.r)
                    !fast index is M for second index of Htmp
                    Htmp=0.0d0
                    atind_mag(1)=lat%cell%ind_mag(neigh%at_pair(1,i_pair))
                    atind_mag(2)=lat%cell%ind_mag(neigh%at_pair(2,i_pair))
                    offset_mag=(atind_mag-1)*3
                    !positive contributions
                    Htmp(offset_mag(1)+2,ind(dim_modes_r,[offset_mag(2)+1,1]))= Hmag*diff_pos(2)     ! Mj_y Mi_x E_x r_y
                    Htmp(offset_mag(1)+3,ind(dim_modes_r,[offset_mag(2)+1,1]))= Hmag*diff_pos(3)     ! Mj_z Mi_x E_x r_z
                    Htmp(offset_mag(1)+1,ind(dim_modes_r,[offset_mag(2)+2,2]))= Hmag*diff_pos(1)     ! Mj_x Mi_y E_y r_x
                    Htmp(offset_mag(1)+3,ind(dim_modes_r,[offset_mag(2)+2,2]))= Hmag*diff_pos(3)     ! Mj_z Mi_y E_y r_z
                    Htmp(offset_mag(1)+1,ind(dim_modes_r,[offset_mag(2)+3,3]))= Hmag*diff_pos(1)     ! Mj_x Mi_z E_z r_x
                    Htmp(offset_mag(1)+2,ind(dim_modes_r,[offset_mag(2)+3,3]))= Hmag*diff_pos(2)     ! Mj_y Mi_z E_z r_y
                    !negative contributions
                    Htmp(offset_mag(1)+1,ind(dim_modes_r,[offset_mag(2)+2,1]))=-Hmag*diff_pos(2)     ! Mj_x Mi_y E_x r_y
                    Htmp(offset_mag(1)+1,ind(dim_modes_r,[offset_mag(2)+3,1]))=-Hmag*diff_pos(3)     ! Mj_x Mi_z E_x r_z
                    Htmp(offset_mag(1)+2,ind(dim_modes_r,[offset_mag(2)+1,2]))=-Hmag*diff_pos(1)     ! Mj_y Mi_x E_y r_x
                    Htmp(offset_mag(1)+2,ind(dim_modes_r,[offset_mag(2)+3,2]))=-Hmag*diff_pos(3)     ! Mj_y Mi_z E_y r_z
                    Htmp(offset_mag(1)+3,ind(dim_modes_r,[offset_mag(2)+1,3]))=-Hmag*diff_pos(1)     ! Mj_z Mi_x E_z r_x
                    Htmp(offset_mag(1)+3,ind(dim_modes_r,[offset_mag(2)+2,3]))=-Hmag*diff_pos(2)     ! Mj_z Mi_y E_z r_y

                    Htmp=-Htmp !flip sign corresponding to previous implementation
                    Call get_coo(Htmp,val_tmp,ind_tmp)

                    !fill Hamiltonian type
                    connect_bnd(2)=neigh%ishell(i_pair)
                    Call Ham_tmp%init_mult_connect_2(neigh%pairs(:,connect_bnd(1):connect_bnd(2)),val_tmp,ind_tmp,"M","ME",lat,2)
                    deallocate(val_tmp,ind_tmp)
                    Call Ham%add(Ham_tmp)
                    Call Ham_tmp%destroy()
                    connect_bnd(1)=connect_bnd(2)+1
                enddo 
            enddo
        enddo
        Ham%desc="antisymmetric magnetoelectric coupling"
        Call mode_set_rank1(Ham%mode_l,"M")
        Call mode_set_rankN(Ham%mode_r,"ME",lat,1)
    endif
end subroutine 
end module 
