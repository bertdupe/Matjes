module m_3spin
use, intrinsic :: iso_fortran_env, only : output_unit,error_unit
use m_input_H_types, only: io_H_sp3
use m_derived_types, only: lattice
use m_neighbor_pair, only: pair_dat_t, get_pair_dat_M

implicit none
character(len=*),parameter  :: ham_desc="3-spin interaction"
private
public :: read_sp3_input, get_3spin

contains

subroutine read_sp3_input(io_param,fname,io)
    use m_io_utils
    integer,intent(in)              :: io_param
    character(len=*), intent(in)    :: fname
    type(io_H_sp3),intent(out)      :: io
    character(len=*),parameter      :: var_name="M_3spin"

    Call get_parameter(io_param,fname,'M_3spin',io%triplet,io%is_set)
    Call get_parameter(io_param,fname,'c_H_sp3',io%c_H_sp3)

end subroutine

subroutine get_3spin(Ham,io,lat)
    !get 4-spin interaction in t_H Hamiltonian format
    use m_H_public, only: t_H, get_Htype
    use m_mode_public
    use m_coo_mat
    use m_neighbor_type, only: neighbors
    use m_setH_util,only: get_coo,ind

    class(t_H),intent(inout)    :: Ham
    type(io_H_sp3),intent(in)   :: io
    type(lattice),intent(in)    :: lat

    !local Hamiltonian
    real(8),allocatable  :: Htmp(:,:)   !local Hamiltonian in (dimmode(1),dimmode(2))-basis
    !local Hamiltonian in coo format
    real(8),allocatable  :: val_tmp(:)
    integer,allocatable  :: ind_tmp(:,:)

    class(t_H),allocatable    :: Ham_tmp    !temporary Hamiltonian type used to add up Ham

    integer         :: i_attrip,N_attrip    !loop parameters which atom-type connection are considered (different neighbor types)
    integer         :: i_dist,N_dist        !loop parameters which  connection are considered (different neighbor types)
    integer         :: i_trip           !loop keeping track which unique connection between the same atom types is considered (indexes "number shells" in neighbors-type)
    integer         :: i_shell          !counting the number of unique connection for given atom types and a distance
    integer         :: connect_bnd(2)   !indices keeping track of which pairs are used for the particular connection
    type(neighbors) :: neigh            !all neighbor information for a given atom-type trip
    real(8)         :: Hmag             !magnitude of Hamiltonian parameter
    integer         :: atind_mag(2)     !index of considered atom in basis of magnetic atoms (1:Nmag)
    integer         :: offset_mag(2)    !offset for start in dim_mode of chosed magnetic atom
    integer         :: dim_modes_r(2)   !mode dimension on the right side (M,M)

    type(coo_mat)   :: mat(2)       !mode construction matrices for rank2 part of Hamiltonian

    if(io%is_set)then
        write(output_unit,'(/2A)') "Start setting Hamiltonian: ", ham_desc
        Call get_Htype(Ham_tmp)
        N_attrip=size(io%triplet)
        allocate(Htmp(lat%M%dim_mode,lat%M%dim_mode**2))!local Hamiltonian modified for each shell/neighbor
        dim_modes_r=[lat%M%dim_mode,lat%M%dim_mode]

        do i_attrip=1,N_attrip
            !loop over different connected atom types
            Call neigh%get(io%triplet(i_attrip)%attype,io%triplet(i_attrip)%dist,lat)
            !write information out
            Call io%triplet(i_attrip)%prt(output_unit,'2X')
            Call neigh%prt(output_unit,'2X')
            stop
            N_dist=size(io%triplet(i_attrip)%dist)
            i_trip=0
            connect_bnd=1 !initialization for lower bound
            do i_dist=1,N_dist
                !loop over distances (nearest, next nearest,... neighbor)
                Hmag=io%triplet(i_attrip)%val(i_dist)*io%c_H_sp3
                do i_shell=1,neigh%Nshell(i_dist)
                    !loop over all different connections with the same distance
                    i_trip=i_trip+1
                    connect_bnd(2)=neigh%ishell(i_trip)

                    !set local Hamiltonian in basis of magnetic orderparameter
                    atind_mag(1)=lat%cell%ind_mag(neigh%at_pair(1,i_trip))
                    atind_mag(2)=lat%cell%ind_mag(neigh%at_pair(2,i_trip))
                    offset_mag=(atind_mag-1)*3
                    Htmp=0.0d0

                    Htmp(offset_mag(1)+1,offset_mag(2)+2)= Hmag
                    Htmp(offset_mag(1)+3,offset_mag(2)+1)= Hmag
                    Htmp(offset_mag(1)+2,offset_mag(2)+3)= Hmag

                    Htmp(offset_mag(1)+2,offset_mag(2)+1)= Hmag
                    Htmp(offset_mag(1)+1,offset_mag(2)+3)= Hmag
                    Htmp(offset_mag(1)+3,offset_mag(2)+2)= Hmag

                    Call get_coo(Htmp,val_tmp,ind_tmp)


                   !fill Hamiltonian type
                    Call Ham_tmp%init_mult_connect_2(neigh%pairs(:,connect_bnd(1):connect_bnd(2)),val_tmp,ind_tmp,"M","MM",lat,3)
                    deallocate(val_tmp,ind_tmp)
                    Call Ham%add(Ham_tmp)
                    Call Ham_tmp%destroy()
                    connect_bnd(1)=connect_bnd(2)+1

                enddo
            enddo
        enddo

        Ham%desc=ham_desc
        Call mode_set_rank1(Ham%mode_l,lat,"M")
#if 0
!this is an obsolete implementation to describe the rank2 mode
!        Call mode_set_rankN(Ham%mode_r,"ME",lat,1)
#else
        Call coo_full_unfold(2,lat%Ncell,dim_modes_r,mat)
        Call mode_set_rankN_sparse(Ham%mode_r,"MM",lat,mat,1)
#endif
    endif
stop
end subroutine

end module
