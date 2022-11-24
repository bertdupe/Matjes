module m_neighbor_pair
!module to construct the basic pair_dat_t, which contains information about how the magnetic nearest neighbor sites are connected
!this is needed for the construction of the mode in case of more than 2 sites as in the 4-spin interaction
use m_derived_types, only: lattice
private
public pair_dat_t, get_pair_dat_U, get_pair_dat_M

type    :: pair_dat_t
    integer,allocatable     :: neighbors(:)     !which neighbors (1st,2nd,3rd,..) are considered
    integer                 :: Nat=0            !number of different atoms
    integer                 :: Npair=0          !number of found pairs combinations
    integer,allocatable     :: Nshell(:)        !(size(neighbors))
    integer,allocatable     :: ind_mag(:)       !(Nat) indices of atoms in magnetic basis
    integer,allocatable     :: ind_ph(:)       !(Nat) indices of atoms in phonon basis
    integer,allocatable     :: pairs(:,:)       !(5,Npair) if about pairs (atid1,atid2,dx,dy,dz)
    integer,allocatable     :: Npair_at(:)         !(Nat*size(neighbors) number of atoms per atom unfolded with neighbors
    real(8),allocatable     :: diff_vec(:,:)    !(3,Npair) real space difference vector from atom 1 to atom2 in real-space 
end type

contains

subroutine get_pair_dat_M(lat,at_type,neigh_in,pair_dat)
    use m_neighbor_type, only: get_neigh_distances
    use, intrinsic :: iso_fortran_env, only : error_unit
    type(lattice),intent(in)                    :: lat
    integer,intent(in)                          :: neigh_in(:)  !wanted neighbors (1=nearest, 2=next nearest,etc...)
    integer,intent(in)                          :: at_type(:)   !atom type indices for which the same-type interaction is to be found
    type(pair_dat_t),allocatable,intent(inout)  :: pair_dat(:)  !output information

    real(8),allocatable     :: atpos1(:,:),atpos2(:,:)

    integer,allocatable     :: Nat_type(:)     !number of atoms per considered type
    integer,allocatable     :: atid(:)         !atom indices of all atoms considered
    integer,allocatable     :: atid_local(:)   !local atom indices 

    integer     :: i_neigh, N_neigh
    integer     :: neigh(size(neigh_in))
    integer     :: iat,Nat    
    integer     :: i_attype, N_attype
    integer     :: i
    integer     :: bnd(2)

    !neighbor parameters
    real(8),allocatable     :: distance(:)
    integer,allocatable     :: pairs(:,:)
    integer,allocatable     :: Nshell(:)
    real(8),allocatable     :: diff_vec(:,:) !difference vector between each shell (3,1:number shells) 
    
    N_attype=size(at_type)
    N_neigh=size(neigh_in)
    neigh=neigh_in+1    !plus 1 to take care for onsite term
    allocate(distance(size(neigh)),source=0.0d0)

    !initial checks
    if(any(at_type>lat%cell%N_attype))then
        write(error_unit,'(3/A,I6,A,I6,A)') "Trying to get neighbor pair for atom type ",maxval(at_type),', but only ',lat%cell%N_attype,' atom-types exist in the cell'
        write(error_unit,'(A)') "Check Hamiltonian input"
        ERROR STOP "ABORT"
    endif
    do i_attype=1,size(at_type) 
        if(any(at_type(i_attype)==at_type(:i_attype-1)).or.any(at_type(i_attype)==at_type(i_attype+1:))) ERROR STOP "EVERY ATOM-TYPE SHALL ONLY APPEAR ONCE in the 4-spin input"
    enddo
    !get atom indices for each atom type
    allocate(Nat_type(N_attype))
    Nat_type=[(count(lat%cell%atomic(:)%type_id==at_type(i_attype)),i_attype=1,N_attype)]
    allocate(atid(sum(Nat_type)),source=0)
    do i_attype=1,N_attype
        bnd=[sum(Nat_type(:i_attype-1))+1,sum(Nat_type(:i_attype))]
        atid(bnd(1):bnd(2))=pack([(i,i=1,size(lat%cell%atomic))],mask=lat%cell%atomic(:)%type_id==at_type(i_attype))
    enddo
    !get actual data
    allocate(pair_dat(N_attype))
    do i_attype=1,N_attype
        Nat=Nat_type(i_attype)

        allocate(atid_local,source=atid(sum(Nat_type(:i_attype-1))+1:sum(Nat_type(:i_attype))))
        allocate(atpos1(3,Nat),atpos2(3,Nat),source=0.0d0)

        do iat=1,Nat
            atpos1(:,iat)=lat%cell%atomic(atid_local(iat))%position
        enddo
        atpos2=atpos1
        !get all neighbor connections(pair_ind,Nshell) and the distances
        Call get_neigh_distances(atpos1, atpos2, neigh, lat, pairs, Nshell, distance, diff_vec) ![2] because onsite would be [1]

        !sort the pairs to be in order of atom types and diff_vec angle
        allocate(pair_dat(i_attype)%Npair_at(N_neigh*Nat),source=0)
        do i_neigh=1,N_neigh
            bnd=[sum(Nshell(:i_neigh-1))+1,sum(Nshell(:i_neigh))]
            Call sort_pairs(Nshell(i_neigh),Nat,pairs(:,bnd(1):bnd(2)),diff_vec(:,bnd(1):bnd(2)),pair_dat(i_attype)%Npair_at((i_neigh-1)*Nat+1:i_neigh*Nat))
        enddo

        !save data in respective pair_dat
        pair_dat(i_attype)%Nat=Nat
        pair_dat(i_attype)%Npair=size(pairs,2)
        allocate(pair_dat(i_attype)%ind_mag,source=atid_local)
        do i=1,size(atid_local)
            pair_dat(i_attype)%ind_mag(i)=lat%cell%ind_mag(pair_dat(i_attype)%ind_mag(i))
        enddo
        Call move_alloc(pairs,pair_dat(i_attype)%pairs)
        Call move_alloc(diff_vec,pair_dat(i_attype)%diff_vec)
        Call move_alloc(Nshell,pair_dat(i_attype)%Nshell)
        deallocate(atpos1,atpos2,atid_local)
    enddo
end subroutine

subroutine get_pair_dat_U(lat,at_type,neigh_in,pair_dat)
    use m_neighbor_type, only: get_neigh_distances
    use, intrinsic :: iso_fortran_env, only : error_unit
    type(lattice),intent(in)                    :: lat
    integer,intent(in)                          :: neigh_in(:)  !wanted neighbors (1=nearest, 2=next nearest,etc...)
    integer,intent(in)                          :: at_type(:)   !atom type indices for which the same-type interaction is to be found
    type(pair_dat_t),allocatable,intent(inout)  :: pair_dat(:)  !output information

    real(8),allocatable     :: atpos1(:,:),atpos2(:,:)

    integer,allocatable     :: Nat_type(:)     !number of atoms per considered type
    integer,allocatable     :: atid(:)         !atom indices of all atoms considered
    integer,allocatable     :: atid_local(:)   !local atom indices

    integer     :: i_neigh, N_neigh
    integer     :: neigh(size(neigh_in))
    integer     :: iat,Nat
    integer     :: i_attype, N_attype
    integer     :: i
    integer     :: bnd(2)

    !neighbor parameters
    real(8),allocatable     :: distance(:)
    integer,allocatable     :: pairs(:,:)
    integer,allocatable     :: Nshell(:)
    real(8),allocatable     :: diff_vec(:,:) !difference vector between each shell (3,1:number shells)

    N_attype=size(at_type)
    N_neigh=size(neigh_in)
    neigh=neigh_in+1    !plus 1 to take care for onsite term
    allocate(distance(size(neigh)),source=0.0d0)

    !initial checks
    if(any(at_type>lat%cell%N_attype))then
        write(error_unit,'(3/A,I6,A,I6,A)') "Trying to get neighbor pair for atom type ",maxval(at_type),', but only ',lat%cell%N_attype,' atom-types exist in the cell'
        write(error_unit,'(A)') "Check Hamiltonian input"
        ERROR STOP "ABORT"
    endif
    do i_attype=1,size(at_type)
        if(any(at_type(i_attype)==at_type(:i_attype-1)).or.any(at_type(i_attype)==at_type(i_attype+1:))) ERROR STOP "EVERY ATOM-TYPE SHALL ONLY APPEAR ONCE in the 4-spin input"
    enddo
    !get atom indices for each atom type
    allocate(Nat_type(N_attype))
    Nat_type=[(count(lat%cell%atomic(:)%type_id==at_type(i_attype)),i_attype=1,N_attype)]
    allocate(atid(sum(Nat_type)),source=0)
    do i_attype=1,N_attype
        bnd=[sum(Nat_type(:i_attype-1))+1,sum(Nat_type(:i_attype))]
        atid(bnd(1):bnd(2))=pack([(i,i=1,size(lat%cell%atomic))],mask=lat%cell%atomic(:)%type_id==at_type(i_attype))
    enddo
    !get actual data
    allocate(pair_dat(N_attype))
    do i_attype=1,N_attype
        Nat=Nat_type(i_attype)
        allocate(atid_local,source=atid(sum(Nat_type(:i_attype-1))+1:sum(Nat_type(:i_attype))))
        allocate(atpos1(3,Nat),atpos2(3,Nat),source=0.0d0)
        do iat=1,Nat
            atpos1(:,iat)=lat%cell%atomic(atid_local(iat))%position
        enddo
        atpos2=atpos1
        !get all neighbor connections(pair_ind,Nshell) and the distances
        Call get_neigh_distances(atpos1, atpos2, neigh, lat, pairs, Nshell, distance, diff_vec) ![2] because onsite would be [1]
        !sort the pairs to be in order of atom types and diff_vec angle
        allocate(pair_dat(i_attype)%Npair_at(N_neigh*Nat),source=0)
        do i_neigh=1,N_neigh
            bnd=[sum(Nshell(:i_neigh-1))+1,sum(Nshell(:i_neigh))]
            Call sort_pairs(Nshell(i_neigh),Nat,pairs(:,bnd(1):bnd(2)),diff_vec(:,bnd(1):bnd(2)),pair_dat(i_attype)%Npair_at((i_neigh-1)*Nat+1:i_neigh*Nat))
        enddo

        !save data in respective pair_dat
        pair_dat(i_attype)%Nat=Nat
        pair_dat(i_attype)%Npair=size(pairs,2)
        allocate(pair_dat(i_attype)%ind_ph,source=atid_local)
        do i=1,size(atid_local)
            pair_dat(i_attype)%ind_ph(i)=lat%cell%ind_ph(pair_dat(i_attype)%ind_ph(i))
        enddo
        Call move_alloc(pairs,pair_dat(i_attype)%pairs)
        Call move_alloc(diff_vec,pair_dat(i_attype)%diff_vec)
        Call move_alloc(Nshell,pair_dat(i_attype)%Nshell)
        deallocate(atpos1,atpos2,atid_local)
    enddo
end subroutine

subroutine sort_pairs(N,Nat,pairs,diff_vec,Npair_at)
    !sort the pairs by the first atom index and then by the in-plane angle 
    use m_sort,only : sort
    use m_constants, only : pi
    integer,intent(in)                  :: N    !number of pairs
    integer,intent(in)                  :: Nat  !number of unique atoms
    integer,intent(inout)               :: pairs(5,N)
    real(8),intent(inout)               :: diff_vec(3,N)
    integer,intent(out)                 :: Npair_at(Nat)

    integer     :: ind(N)
    integer     :: tmp(N)
    integer     :: iat
    integer     :: bnd(2)
    real(8)     :: tmp_vec(3,N)
    real(8)     :: angle(N)

    do iat=1,Nat
       Npair_at(iat)=count(pairs(1,:)==iat)
    enddo
    tmp=pairs(1,:)
    ind=0
    Call sort(N,tmp,ind)    !sort by first atom index
    do iat=1,Nat
        bnd=[sum(Npair_at(:iat-1))+1,sum(Npair_at(:iat))]
        tmp_vec(:,1:Npair_at(iat))=diff_vec(:,ind(bnd(1):bnd(2)))
        angle(1:Npair_at(iat))=atan2(tmp_vec(2,1:Npair_at(iat)),tmp_vec(1,1:Npair_at(iat)))
        Call sort(Npair_at(iat),angle,ind(bnd(1):bnd(2)),1.0d-6)    !sort by angle where first atom index is the same
    enddo
    !apply permutations
    pairs=pairs(:,ind)
    diff_vec=diff_vec(:,ind)
end subroutine
end module

