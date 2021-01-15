module m_exchange_heisenberg_D
use m_input_H_types, only: io_H_D
implicit none
private
public read_D_input, get_exchange_D
contains
subroutine read_D_input(io_param,fname,io)
    use m_io_utils
    integer,intent(in)              :: io_param
    character(len=*), intent(in)    :: fname
    type(io_H_D),intent(out)        :: io

    Call get_parameter(io_param,fname,'magnetic_D',io%trip,io%is_set) 
end subroutine


subroutine get_exchange_D(Ham,io,lat)
    !get coupling in t_H Hamiltonian format
    use m_H_public
    use m_derived_types
    use m_setH_util,only: get_coo
    use m_neighbor_type, only: neighbors

    class(t_H),intent(inout)    :: Ham
    type(io_H_D),intent(in)     :: io
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
    real(8)         :: DMI(3)           !local DMI vector
    logical         :: is_set           

    if(io%is_set)then
        Call get_Htype(Ham_tmp)
        N_attrip=size(io%trip)
        allocate(Htmp(lat%M%dim_mode,lat%M%dim_mode))!local Hamiltonian modified for each shell/neighbor
        do i_attrip=1,N_attrip
            !loop over different connected atom types
            Call neigh%get(io%trip(i_attrip)%attype(1:2),io%trip(i_attrip)%dist,lat)
            N_dist=size(io%trip(i_attrip)%dist)
            i_trip=0
            connect_bnd=1 !initialization for lower bound
            do i_dist=1,N_dist
                !loop over distances (nearest, next nearest,... neighbor)
                Hmag=io%trip(i_attrip)%val(i_dist)
                do i_shell=1,neigh%Nshell(i_dist)
                    !loop over all different connections with the same distance
                    i_trip=i_trip+1
                    Call get_DMI(neigh%at_pair(:,i_trip),neigh%pairs(:,connect_bnd(1)),io%trip(i_attrip)%attype(3),lat,DMI)
                    Call print_DMI(DMI,neigh%at_pair(:,i_trip),neigh%pairs(:,connect_bnd(1)),io%trip(i_attrip)%attype(3),io%trip(i_attrip)%dist(i_dist),lat)
                    connect_bnd(2)=neigh%ishell(i_trip)
                    if(norm2(DMI)>1.0d-8)then 
                        !set local Hamiltonian in basis of magnetic orderparameter
                        atind_mag(1)=lat%cell%ind_mag(neigh%at_pair(1,i_trip))
                        atind_mag(2)=lat%cell%ind_mag(neigh%at_pair(2,i_trip))
                        offset_mag=(atind_mag-1)*3
                        Htmp=0.0d0

                        Htmp(offset_mag(1)+1,offset_mag(2)+2)= DMI(3) * Hmag
                        Htmp(offset_mag(1)+3,offset_mag(2)+1)= DMI(2) * Hmag
                        Htmp(offset_mag(1)+2,offset_mag(2)+3)= DMI(1) * Hmag

                        Htmp(offset_mag(1)+2,offset_mag(2)+1)=-DMI(3) * Hmag
                        Htmp(offset_mag(1)+1,offset_mag(2)+3)=-DMI(2) * Hmag
                        Htmp(offset_mag(1)+3,offset_mag(2)+2)=-DMI(1) * Hmag
                        Htmp=-Htmp !flip sign corresponding to previous implementation

                        Call get_coo(Htmp,val_tmp,ind_tmp)

                        !fill Hamiltonian type
                        Call Ham_tmp%init_connect(neigh%pairs(:,connect_bnd(1):connect_bnd(2)),val_tmp,ind_tmp,"MM",lat,2)
                        deallocate(val_tmp,ind_tmp)
                        Call Ham%add(Ham_tmp)
                        Call Ham_tmp%destroy()
                    endif
                    connect_bnd(1)=connect_bnd(2)+1
                enddo 
            enddo
        enddo

        Ham%desc="antisymmetric magnetic exchange"
        if(.not.Ham%is_set())then
            write(*,'(/A)') "DID NOT FIND A SINGLE DMI CONTRIBUTION, LEAVE THE DMI OUT IF YOU EXPECT THE DMI TO VANISH OR FIND A MISTAKE IN THIS CALCULATION"
            !continuing will give problems since the t_H-type Ham is not set but most probably will be assumed to be set in routine calling this subroutine
            STOP "CHECK INPUT"
        endif
    endif
end subroutine 


subroutine get_DMI(atom_mag,pair_mag,atom_get_type,lat,DMI_sum)
    !Temporary way to get basic working version with DMI in case of new input with several magnetic/non-magnetic atoms.
    !The idea is that you give 2 magnetic atoms specified by atom_mag and pair_mag and the type of the non-magnetic atom causing the DMI (atom_get_type).
    !Then atoms of the DMI-causing atomtype that are nearest to the center of the magnetic atoms are searched.
    !Those atoms are assumed to cause DMI with normalized vectors, which is summed into DMI_sum and normalized again.

    !WOULD IT NOT RATHER MAKE SENSE TO CALCULATE THE DM-CAUSING ATOMS RELATIVE TO THE MOST CENTERED UNITCELL TO AVOID PERIODICITY PROBLEMS?
    use m_derived_types, only: lattice
    use m_neighbor_type, only: get_neigh_distances
    use m_constants, only : pi
    use m_vector, only: cross
    integer,intent(in)          :: atom_mag(2)      !both connected magnetic atoms in (1:Nat)-space
    integer,intent(in)          :: pair_mag(2)      !indices of connected magnetic atoms in (1:Ncell)-space
    integer,intent(in)          :: atom_get_type    !atom type of the non-magnetic source atom
    type(lattice),intent(in)    :: lat              !allmighty lattice type with all knowledge over the geometry
    real(8),intent(out)         :: DMI_sum(3)       !sum of all DMI-vectors 

    integer     :: ind4_mag(4,2)    !indices of both magnetic atoms is (ix,iy,iz,ia)-space
    real(8)     :: pos_mag(3,2)     !position of magnetic atoms in real-space
    real(8)     :: pos_center(3)    !center between magnetic atoms in real-space
    real(8)     :: pos_center_uc(3) !center in the first unit_cell
    real(8)     :: diff(3)          !temporary difference vector
    integer     :: ind_cell(3)      !indices for moving pos_center_uc to the first unit_cell

    integer,allocatable     :: id_nonM(:)       !indices of non-magnetic atoms
    real(8),allocatable     :: atpos_nonM(:,:)  !position in atoms in first unit-cell of atom_get_type type
    integer,allocatable     :: pairs(:,:)       !
    integer,allocatable     :: Nshell(:)        !(1)-array which gives the number of atoms at same distance from center-position
    real(8)                 :: distance(1)      !distance from magnetic atom
    real(8)                 :: pos_found(3)     !position of non-magnetic atom from which a DMI is calculated
    real(8)                 :: pos_diff(3,2)    !difference vec between magnetic atom positions and found non-magnetic atom

    real(8)                 :: DMI(3)   !DMI as caused by each found non-magnetic atom
    real(8)                 :: nrm      !norm

    integer     ::  i

    !get real-space positions of magnetic atoms and center(obeying supercell-symmetry)
    ind4_mag(1:3,1)=lat%index_1_3(pair_mag(1)) 
    ind4_mag(1:3,2)=lat%index_1_3(pair_mag(2)) 
    ind4_mag(4,:)=atom_mag
    Call lat%pos_ind(ind4_mag(:,1),pos_mag(:,1))
    Call lat%pos_ind(ind4_mag(:,2),pos_mag(:,2))
    diff=pos_mag(:,2)-pos_mag(:,1)
    Call lat%min_diffvec(diff)
    pos_mag(:,2)=diff+pos_mag(:,1)
    pos_center=sum(pos_mag,2)*0.5d0
    !get pos_center_uc 
    ind_cell=floor(matmul(pos_center,lat%astar)/(2.0d0*pi))
    pos_center_uc=pos_center-matmul(real(ind_cell,8),lat%areal)

    !fill at_pos_nonM with all atoms of the type in uc
    Call lat%cell%ind_attype(atom_get_type,id_nonM)
    allocate(atpos_nonM(3,size(id_nonM)))
    do i=1,size(id_nonM)
        atpos_nonM(:,i)=lat%cell%atomic(id_nonM(i))%position
    enddo
    !get all considered non-magnetic atoms with minimal distance from pos_center_uc
    Call get_neigh_distances(atpos_nonM,reshape(pos_center_uc,[3,1]),[1],lat,pairs,Nshell,distance)

    !sum the DMI caused by each found non-magnetic atom
    DMI_sum=0.0d0
    do i=1,Nshell(1)
        pos_found=lat%cell%atomic(id_nonM(pairs(1,i)))%position
        pos_found=pos_found+matmul(real(pairs(3:5,i),8),lat%areal)
        pos_found=pos_found-pos_center_uc+pos_center  !finally position of considered non-magnetic atom

        pos_diff(:,1)=pos_mag(:,1)-pos_found
        pos_diff(:,2)=pos_mag(:,2)-pos_found
        DMI=cross(pos_diff(:,1),pos_diff(:,2))
        nrm=norm2(DMI)
        if(nrm>1.0d-8) DMI_sum=DMI_sum+DMI/nrm
    enddo
    nrm=norm2(DMI_sum)
    if(nrm>1.0d-8) DMI_sum=DMI_sum/nrm
end subroutine

subroutine print_DMI(DMI,at_pair,pairs,type_mediate,dist,lat)
    use, intrinsic :: iso_fortran_env,only : output_unit
    use m_derived_types, only: lattice
    real(8),intent(in)          :: DMI(3)
    integer,intent(in)          :: at_pair(2)
    integer,intent(in)          :: pairs(2) 
    integer,intent(in)          :: type_mediate
    integer,intent(in)          :: dist
    type(lattice),intent(in)    :: lat

    integer         ::  transl(3)
    integer         ::  i

    transl=lat%index_1_3(pairs(2))-lat%index_1_3(pairs(1))
    !make sure periodicity is used for shortes translation
    do i=1,3
        if(lat%periodic(i))then
            if(abs(transl(i)-lat%dim_lat(i))<abs(transl(i))) transl(i)=transl(i)-lat%dim_lat(i)
            if(abs(transl(i)+lat%dim_lat(i))<abs(transl(i))) transl(i)=transl(i)+lat%dim_lat(i)
        endif
    enddo

    write(output_unit,'(/A,3(/F12.6))') 'Found DMI-vector',DMI
    write(output_unit,'(A,2I6)')        '  between atoms:',at_pair
    write(output_unit,'(A,3I6)')        '  with lattice translation:',transl
    write(output_unit,'(A,I6)')         '  connected by atom type:',type_mediate
    write(output_unit,'(A,I4,A/)')      '  corresponding to distance no.', dist
end subroutine
end module
