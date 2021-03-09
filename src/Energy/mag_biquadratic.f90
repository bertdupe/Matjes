module m_Mag_Biq
use m_input_H_types, only: io_H_Mag_Biq
use m_io_utils, only: get_parameter,get_coeff,number_nonzero_coeff,max_ind_variable
implicit none
private
public :: get_Mag_Biq, read_Mag_Biq_input
contains

subroutine read_Mag_Biq_input(io_param,fname,io)
    integer,intent(in)              :: io_param
    character(len=*), intent(in)    :: fname
    type(io_H_Mag_Biq),intent(out)  :: io

    Call get_parameter(io_param,fname,'M_biq',io%pair,io%is_set)
end subroutine

subroutine get_Mag_Biq(Ham,io,lat)
    !get coupling in t_H Hamiltonian format
    use m_H_public
    use m_derived_types, only: lattice
    use m_setH_util, only: get_coo
    use m_neighbor_type, only: neighbors

    class(t_H),intent(inout)        :: Ham   !Hamiltonian in which all contributions are added up
    type(io_H_Mag_Biq),intent(in)   :: io
    type(lattice),intent(in)        :: lat

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
    real(8)         :: Biq              !magnitude of Hamiltonian parameter
    integer         :: atind_mag(2)     !index of considered atom in basis of magnetic atoms (1:Nmag)

    if(io%is_set)then
        Call get_Htype(Ham_tmp)
        N_atpair=size(io%pair)
        allocate(Htmp(lat%M%dim_mode,lat%M%dim_mode))!local Hamiltonian modified for each shell/neighbor
        do i_atpair=1,N_atpair
            !loop over different connected atom types
            Call neigh%get(io%pair(i_atpair)%attype,io%pair(i_atpair)%dist,lat)
            N_dist=size(io%pair(i_atpair)%dist)
            i_pair=0
            connect_bnd=1 !initialization for lower bound
            do i_dist=1,N_dist
                !loop over distances (nearest, next nearest,... neighbor)
                Biq=io%pair(i_atpair)%val(i_dist)
                do i_shell=1,neigh%Nshell(i_dist)
                    !loop over all different connections with the same distance
                    i_pair=i_pair+1

                    !set local Hamiltonian in basis of magnetic orderparameter
                    atind_mag(1)=lat%cell%ind_mag(neigh%at_pair(1,i_pair))
                    atind_mag(2)=lat%cell%ind_mag(neigh%at_pair(2,i_pair))
                    Htmp=0.0d0
                    Htmp(atind_mag(1)*3-2,atind_mag(2)*3-2)=Biq
                    Htmp(atind_mag(1)*3-1,atind_mag(2)*3-1)=Biq
                    Htmp(atind_mag(1)*3  ,atind_mag(2)*3  )=Biq
                    connect_bnd(2)=neigh%ishell(i_pair)
                    Htmp=-Htmp !flip sign corresponding to previous implementation
                    Call get_coo(Htmp,val_tmp,ind_tmp)

                    !fill Hamiltonian type
                    Call Ham_tmp%init_mult_connect_2(neigh%pairs(:,connect_bnd(1):connect_bnd(2)),val_tmp,ind_tmp,"MM","MM",lat,2)
                    deallocate(val_tmp,ind_tmp)
                    Call Ham%add(Ham_tmp)
                    Call Ham_tmp%destroy()
                    connect_bnd(1)=connect_bnd(2)+1
                enddo
            enddo
        enddo
        Ham%desc="Biquadratic magnetic exchange"
    endif
end subroutine

end module m_Mag_Biq
