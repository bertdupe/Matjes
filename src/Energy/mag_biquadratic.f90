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


#if 0
subroutine get_Mag_Biq(Ham,io,lat)
    !get coupling in t_H Hamiltonian format
    use m_H_public
    use m_derived_types, only: lattice
    use m_setH_util, only: get_coo
    use m_neighbor_type, only: neighbors
    use m_mode_public

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
    integer         :: ind_ham(3,2)

    if(io%is_set)then
        Call get_Htype(Ham_tmp)
        N_atpair=size(io%pair)
        allocate(Htmp(lat%M%dim_mode**2,lat%M%dim_mode**2))!local Hamiltonian modified for each shell/neighbor
        Nmag=lat%nmag
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
                    ind_ham(:,1)=get_diag_MM(atind_mag(1),Nmag)
                    ind_ham(:,2)=get_diag_MM(atind_mag(2),Nmag)
                    Htmp=0.0d0
                    Htmp(ind_ham(1,1),ind_ham(1,2))=Biq
                    Htmp(ind_ham(2,1),ind_ham(2,2))=Biq
                    Htmp(ind_ham(3,1),ind_ham(3,2))=Biq
                    connect_bnd(2)=neigh%ishell(i_pair)
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
        Call mode_set_rankN(Ham%mode_l,"MM",lat,1)
        Call mode_set_rankN(Ham%mode_r,"MM",lat,1)
    endif
end subroutine

function get_diag_MM(ind_mag,Nmag)result(ind)
    !function which returns the 3 matrix indices for the diagonal MxM entries of a fully unfolded MxX mode 
    !(i.e the M_ind_mag_x * M_ind_mag_x, M_ind_mag_y * M_ind_mag_y, pM_ind_mag_z * M_ind_mag_z
    integer,intent(in)  :: ind_mag  ![1,Nmag], which magnetic atom is consideres
    integer,intent(in)  :: Nmag     !how many magnetic atoms are in one unit-cell
    integer             :: ind(3)

    ind=(ind_mag-1)*Nmag*9  !get to offset where the slower index starts with the correct ind_mag offset
    ind=ind+(ind_mag-1)*3   !advance the faster index to get to the correct index start offset
    ind=ind(1)+[1,3*Nmag+2,6*Nmag+3]    !add offsets to get xx,yy,zz index
end function 


#else
!faster sparse implementation

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
    integer         :: ind_ham(3,2)
    !parameter to construct the modes
    type(coo_mat)   :: mat(2)   !mode construction matrices of left/right side of Hamiltonian ( first left, reused for right)
    integer         :: dim_mat(2),nnz
    integer,allocatable ::  row(:),col(:)
    real(8),allocatable ::  val(:)
    integer         :: i,j,ii

    if(io%is_set)then
        Call get_Htype(Ham_tmp)
        N_atpair=size(io%pair)
        allocate(Htmp(lat%M%dim_mode,lat%M%dim_mode))!local Hamiltonian modified for each shell/neighbor
        Nmag=lat%nmag
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
                    ind_ham(:,1)=[((atind_mag(1)-1)*3+i,i=1,3)]
                    ind_ham(:,2)=[((atind_mag(2)-1)*3+i,i=1,3)]
                    Htmp=0.0d0
                    Htmp(ind_ham(1,1),ind_ham(1,2))=Biq
                    Htmp(ind_ham(2,1),ind_ham(2,2))=Biq
                    Htmp(ind_ham(3,1),ind_ham(3,2))=Biq
                    connect_bnd(2)=neigh%ishell(i_pair)
                    Call get_coo(Htmp,val_tmp,ind_tmp)

                    !fill Hamiltonian type
                    Call Ham_tmp%init_mult_connect_2(neigh%pairs(:,connect_bnd(1):connect_bnd(2)),val_tmp,ind_tmp,"MM","MM",lat,2,[3*Nmag,3*Nmag])
                    deallocate(val_tmp,ind_tmp)
                    Call Ham%add(Ham_tmp)
                    Call Ham_tmp%destroy()
                    connect_bnd(1)=connect_bnd(2)+1
                enddo
            enddo
        enddo
        Ham%desc="Biquadratic magnetic exchange"

        !first, quicker varying mode
        dim_mat=[lat%M%dim_mode*lat%Ncell,lat%M%dim_mode*lat%Ncell]
        nnz=lat%M%dim_mode*lat%Ncell
        allocate(row(nnz),col(nnz),source=0)
        allocate(val(nnz),source=1.0d0)
        ii=0
        do i=0,lat%ncell-1
            do j=0,Nmag-1
                row(i*Nmag*3+j*3+1)=i*Nmag*3+j*3+1
                row(i*Nmag*3+j*3+2)=i*Nmag*3+j*3+2
                row(i*Nmag*3+j*3+3)=i*Nmag*3+j*3+3
                ii=ii+3
            enddo
        enddo
        col=row !only diagonal terms

        !both M-states have same construction procedure
        Call mat(1)%init(dim_mat,nnz,row,col,val)
        Call mat(2)%init(dim_mat,nnz,row,col,val)
        Call mode_set_rankN_sparse(Ham%mode_l,"MM",lat,mat,1)

        !recreate mat since it is destroyed in initialization
        Call mat(1)%init(dim_mat,nnz,row,col,val)
        Call mat(2)%init(dim_mat,nnz,row,col,val)
        Call mode_set_rankN_sparse(Ham%mode_r,"MM",lat,mat,1)
    endif
end subroutine
#endif

end module m_Mag_Biq
