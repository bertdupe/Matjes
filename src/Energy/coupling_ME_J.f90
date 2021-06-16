module m_coupling_ME_J
use m_input_H_types, only: io_H_ME_J
use m_io_utils, only: get_parameter,get_coeff,number_nonzero_coeff,max_ind_variable
use m_coo_mat
implicit none
private
public :: get_coupling_ME_J, read_ME_J_input
contains

subroutine read_ME_J_input(io_param,fname,io)
    integer,intent(in)              :: io_param
    character(len=*), intent(in)    :: fname
    type(io_H_ME_J),intent(out)       :: io

    Call get_parameter(io_param,fname,'ME_sym',io%pair,io%is_set) 
end subroutine

subroutine get_coupling_ME_J(Ham,io,lat)
    !get coupling in t_H Hamiltonian format
    use m_H_public
    use m_derived_types
    use m_setH_util,only: get_coo,ind
    use m_neighbor_type, only: neighbors
    use m_mode_public

    class(t_H),intent(inout)    :: Ham
    type(io_H_ME_J),intent(in)  :: io
    type(lattice),intent(in)    :: lat

    !local Hamiltonian
    real(8),allocatable  :: Htmp(:,:)   !local Hamiltonian in (dimmode(1),dimmode(2))-basis
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

    integer         :: dim_modes_r(2)   !mode dimension on the right side (M,E)

    type(coo_mat)   :: mat(2)       !mode construction matrices for rank2 part of Hamiltonian


    if(io%is_set)then
        if(lat%E%dim_mode==0) STOP "E-field has to be set when using coupling_ME"
        Call get_Htype(Ham_tmp)
        N_atpair=size(io%pair)
        allocate(Htmp(lat%M%dim_mode,lat%M%dim_mode*lat%E%dim_mode))!local Hamiltonian modified for each shell/neighbor
        dim_modes_r=[lat%M%dim_mode,lat%E%dim_mode]
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
                    Htmp(offset_mag(1)+1,ind(dim_modes_r,[offset_mag(2)+1,1]))=J     ! Mi_x Mj_x E_x
                    Htmp(offset_mag(1)+2,ind(dim_modes_r,[offset_mag(2)+2,2]))=J     ! Mi_y Mj_y E_y
                    Htmp(offset_mag(1)+3,ind(dim_modes_r,[offset_mag(2)+3,3]))=J     ! Mi_z Mj_z E_z
                    connect_bnd(2)=neigh%ishell(i_pair)
                    Htmp=-Htmp !flip sign corresponding to previous implementation
                    Call get_coo(Htmp,val_tmp,ind_tmp)

                    !fill Hamiltonian type
                    Call Ham_tmp%init_mult_connect_2(neigh%pairs(:,connect_bnd(1):connect_bnd(2)),val_tmp,ind_tmp,"M","ME",lat,2)
                    deallocate(val_tmp,ind_tmp)
                    Call Ham%add(Ham_tmp)
                    Call Ham_tmp%destroy()
                    connect_bnd(1)=connect_bnd(2)+1
                enddo 
            enddo
        enddo
        Ham%desc="symmetric magnetoelectric coupling"
        Call mode_set_rank1(Ham%mode_l,lat,"M")
#if 0
!this is an obsolete implementation to describe the rank2 mode
!        Call mode_set_rankN(Ham%mode_r,"ME",lat,1)
#else
        Call coo_full_unfold(2,lat%Ncell,dim_modes_r,mat)
        Call mode_set_rankN_sparse(Ham%mode_r,"ME",lat,mat,1)
#endif
    endif
end subroutine 
end module 
