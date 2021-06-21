module m_exchange_TJ
use m_input_H_types, only: io_H_TJ
implicit none
private
public read_TJ_input, get_exchange_TJ
contains
subroutine read_TJ_input(io_param,fname,io)
    use m_io_utils
    integer,intent(in)              :: io_param
    character(len=*), intent(in)    :: fname
    type(io_H_TJ),intent(out)        :: io

    Call get_parameter(io_param,fname,'MT_J',io%pair,io%is_set) 
end subroutine


subroutine get_exchange_TJ(Ham,io,lat)
    !get coupling in t_H Hamiltonian format
    use m_H_public
    use m_derived_types
    use m_setH_util,only: get_coo
    use m_neighbor_type, only: neighbors
    use m_mode_public
    use m_coo_mat

    class(t_H),intent(inout)    :: Ham
    type(io_H_TJ),intent(in)    :: io
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
    real(8)         :: J                !magnitude of Hamiltonian parameter
    integer         :: atind_mag(2)     !index of considered atom in basis of magnetic atoms (1:Nmag)A
    integer         :: offset_mag(2)    !offset for start in dim_mode of chosed magnetic atom
    type(coo_mat)   :: mat(3)       !mode construction matrices for rank3 part of Hamiltonian

    if(io%is_set)then
        if(lat%T%dim_mode==0) STOP "T-field has to be set when using TJ-exchange"
        Call get_Htype(Ham_tmp)
        N_atpair=size(io%pair)
        if(lat%T%dim_mode/=1) ERROR STOP "THIS HAS TO BE UPDATED IF THE TEMPERATURE IS NO LONGER A SIMPLE SCALAR"
        allocate(Htmp(lat%M%dim_mode,lat%M%dim_mode))!local Hamiltonian modified for each shell/neighbor
        do i_atpair=1,N_atpair
            !loop over different connected atom types
            Call neigh%get(io%pair(i_atpair)%attype,io%pair(i_atpair)%dist,lat)
            N_dist=size(io%pair(i_atpair)%dist)
            i_pair=0
            connect_bnd=1 !initialization for lower bound
            do i_dist=1,N_dist
                !loop over distances (nearest, next nearest,... neighbor)
                J=io%pair(i_atpair)%val(i_dist)
                do i_shell=1,neigh%Nshell(i_dist)
                    !loop over all different connections with the same distance
                    i_pair=i_pair+1

                    !set local Hamiltonian in basis of magnetic orderparameter
                    atind_mag(1)=lat%cell%ind_mag(neigh%at_pair(1,i_pair))
                    atind_mag(2)=lat%cell%ind_mag(neigh%at_pair(2,i_pair))
                    offset_mag=(atind_mag-1)*3
                    Htmp=0.0d0
                    Htmp(offset_mag(1)+1,offset_mag(2)+1)=J     ! Mi_x Mj_x Tj Tj
                    Htmp(offset_mag(1)+2,offset_mag(2)+2)=J     ! Mi_y Mj_y Tj Tj
                    Htmp(offset_mag(1)+3,offset_mag(2)+3)=J     ! Mi_z Mj_z Tj Tj
                    Htmp=-Htmp !flip sign corresponding to previous implementation
                    connect_bnd(2)=neigh%ishell(i_pair)
                    Call get_coo(Htmp,val_tmp,ind_tmp)

                    !fill Hamiltonian type
                    Call Ham_tmp%init_mult_connect_2(neigh%pairs(:,connect_bnd(1):connect_bnd(2)),val_tmp,ind_tmp,"M","MTT",lat,2)
                    deallocate(val_tmp,ind_tmp)
                    Call Ham%add(Ham_tmp)
                    Call Ham_tmp%destroy()
                    connect_bnd(1)=connect_bnd(2)+1
                enddo 
            enddo
        enddo
        Ham%desc="T^2 M^2 exchange"
        Call mode_set_rank1(Ham%mode_l,lat,"M")
#if 0
    !obsolete
!        Call mode_set_rankN(Ham%mode_r,"MTT",lat,1)
#else
        Call coo_full_unfold(3,lat%Ncell,[lat%M%dim_mode,lat%T%dim_mode,lat%T%dim_mode],mat)
        Call mode_set_rankN_sparse(Ham%mode_r,"MTT",lat,mat,1)
#endif
    endif
end subroutine 

end module
