module m_Mag_Biq
use, intrinsic :: iso_fortran_env, only : output_unit
use m_input_H_types, only: io_H_Mag_Biq
use m_io_utils, only: get_parameter,get_coeff,number_nonzero_coeff,max_ind_variable
implicit none
character(len=*),parameter  :: ham_desc="Biquadratic magnetic exchange"
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
    use m_mode_public
    use m_coo_mat

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
    integer         :: Nmag             
    integer         :: ind_ham(9,2)
    !parameter to construct the modes
    type(coo_mat)   :: mat(2)   !mode construction matrices of left/right side of Hamiltonian ( first left, second right)
    integer         :: dim_mat(2),nnz
    integer,allocatable ::  row(:),col(:)
    real(8),allocatable ::  val(:)
    integer         :: i,j,ii

    if(io%is_set)then
        Call get_Htype(Ham_tmp)
        N_atpair=size(io%pair)
        Nmag=lat%nmag
        allocate(Htmp(Nmag*9,Nmag*9))!local Hamiltonian modified for each shell/neighbor
        do i_atpair=1,N_atpair
            !loop over different connected atom types
            Call neigh%get(io%pair(i_atpair)%attype,io%pair(i_atpair)%dist,lat)
            !write information out
            Call io%pair(i_atpair)%prt(output_unit,'2X')
            Call neigh%prt(output_unit,'2X')
            N_dist=size(io%pair(i_atpair)%dist)
            i_pair=0
            connect_bnd=1 !initialization for lower bound
            do i_dist=1,N_dist
                !loop over distances (nearest, next nearest,... neighbor)
                Biq=-io%pair(i_atpair)%val(i_dist)
                do i_shell=1,neigh%Nshell(i_dist)
                    !loop over all different connections with the same distance
                    i_pair=i_pair+1

                    !set local Hamiltonian in basis of magnetic orderparameter
                    atind_mag(1)=lat%cell%ind_mag(neigh%at_pair(1,i_pair))
                    atind_mag(2)=lat%cell%ind_mag(neigh%at_pair(2,i_pair))
                    ind_ham(:,1)=[((atind_mag(1)-1)*3+i,i=1,9)]
                    ind_ham(:,2)=[((atind_mag(2)-1)*3+i,i=1,9)]
                    Htmp=0.0d0
                    do i=1,9
                        Htmp(ind_ham(i,1),ind_ham(i,2))=Biq
                    end do
                    connect_bnd(2)=neigh%ishell(i_pair)
                    Call get_coo(Htmp,val_tmp,ind_tmp)

                    !fill Hamiltonian type
                    Call Ham_tmp%init_mult_connect_2(neigh%pairs(:,connect_bnd(1):connect_bnd(2)),val_tmp,ind_tmp,"MM","MM",lat,2,[9*Nmag,9*Nmag])
                    deallocate(val_tmp,ind_tmp)
                    Call Ham%add(Ham_tmp)
                    Call Ham_tmp%destroy()
                    connect_bnd(1)=connect_bnd(2)+1
                enddo
            enddo
        enddo
        Ham%desc=ham_desc

        !first, quicker varying mode
        dim_mat=[Nmag*9*lat%Ncell,lat%M%dim_mode*lat%Ncell]
        nnz=dim_mat(1)
        allocate(row(nnz),col(nnz),source=0)
        allocate(val(nnz),source=1.0d0)
        row=[(i,i=1,nnz)]
        ii=0
        do i=0,lat%ncell-1
            do j=0,Nmag-1
                col(i*Nmag*9+j*3+1)=i*Nmag*3+j*3+1
                col(i*Nmag*9+j*3+2)=i*Nmag*3+j*3+2
                col(i*Nmag*9+j*3+3)=i*Nmag*3+j*3+3
                col(i*Nmag*9+j*3+4)=i*Nmag*3+j*3+1
                col(i*Nmag*9+j*3+5)=i*Nmag*3+j*3+2
                col(i*Nmag*9+j*3+6)=i*Nmag*3+j*3+3
                col(i*Nmag*9+j*3+7)=i*Nmag*3+j*3+1
                col(i*Nmag*9+j*3+8)=i*Nmag*3+j*3+2
                col(i*Nmag*9+j*3+9)=i*Nmag*3+j*3+3
                ii=ii+9
            enddo
        enddo
        Call mat(1)%init(dim_mat,nnz,row,col,val)
        ii=0
        do i=0,lat%ncell-1
            do j=0,Nmag-1
                col(i*Nmag*9+j*3+1)=i*Nmag*3+j*3+1
                col(i*Nmag*9+j*3+2)=i*Nmag*3+j*3+1
                col(i*Nmag*9+j*3+3)=i*Nmag*3+j*3+1
                col(i*Nmag*9+j*3+4)=i*Nmag*3+j*3+2
                col(i*Nmag*9+j*3+5)=i*Nmag*3+j*3+2
                col(i*Nmag*9+j*3+6)=i*Nmag*3+j*3+2
                col(i*Nmag*9+j*3+7)=i*Nmag*3+j*3+3
                col(i*Nmag*9+j*3+8)=i*Nmag*3+j*3+3
                col(i*Nmag*9+j*3+9)=i*Nmag*3+j*3+3
                ii=ii+9
            enddo
        enddo
        Call mat(2)%init(dim_mat,nnz,row,col,val)
        Call mode_set_rankN_sparse(Ham%mode_l,"MM",lat,mat,1)
        Call Ham%mode_l%copy(Ham%mode_r)
    endif
end subroutine

end module m_Mag_Biq
