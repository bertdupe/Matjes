module m_exchange_heisenberg_D
use, intrinsic :: iso_fortran_env, only : output_unit
use m_input_H_types, only: io_H_D
implicit none
private
public read_D_input, get_exchange_D, get_exchange_D_fft

character(len=*),parameter  :: ham_desc="antisymmetric magnetic exchange"

contains
subroutine read_D_input(io_param,fname,io)
    use m_io_utils
    integer,intent(in)              :: io_param
    character(len=*), intent(in)    :: fname
    type(io_H_D),intent(out)        :: io

    Call get_parameter(io_param,fname,'magnetic_D',io%trip,io%is_set) 
    Call get_parameter(io_param,fname,'magnetic_D_fft',io%fft) 
end subroutine


subroutine get_exchange_D(Ham,io,lat,Ham_shell_pos,neighbor_pos_list)
    !get coupling in t_H Hamiltonian format
    use m_H_public
    use m_derived_types
    use m_setH_util,only: get_coo
    use m_neighbor_type, only: neighbors
    use m_mode_public

    class(t_H),intent(inout)                       :: Ham
    type(io_H_D),intent(in)                        :: io
    type(lattice),intent(in)                       :: lat
    real(8),optional,allocatable,intent(inout)     :: neighbor_pos_list(:,:)
    real(8),optional,allocatable,intent(inout)     :: Ham_shell_pos(:,:,:)

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

    if(io%is_set)then
        write(output_unit,'(/2A)') "Start setting Hamiltonian: ", ham_desc
        Call get_Htype(Ham_tmp)
        N_attrip=size(io%trip)
        allocate(Htmp(lat%M%dim_mode,lat%M%dim_mode))!local Hamiltonian modified for each shell/neighbor

        if (present(Ham_shell_pos)) then
          write(output_unit,'(/2A)') "Preparing the Fourier Transform of Hamiltonian: ", ham_desc
          i_trip=0
          do i_attrip=1,N_attrip
             Call neigh%get(io%trip(i_attrip)%attype,io%trip(i_attrip)%dist,lat)
             N_dist=size(io%trip(i_attrip)%dist)
             do i_dist=1,N_dist
                do i_shell=1,neigh%Nshell(i_dist)
                   i_trip=i_trip+1
                enddo
             enddo
          enddo
          allocate(Ham_shell_pos(lat%M%dim_mode,lat%M%dim_mode,i_trip))
          Ham_shell_pos=0.0d0
          allocate(neighbor_pos_list(3,i_trip))
        endif

        do i_attrip=1,N_attrip
            !loop over different connected atom types
            Call neigh%get(io%trip(i_attrip)%attype(1:2),io%trip(i_attrip)%dist,lat)
            !write information out
            Call io%trip(i_attrip)%prt(output_unit,'2X')
            Call neigh%prt(output_unit,'2X')
            N_dist=size(io%trip(i_attrip)%dist)
            i_trip=0
            connect_bnd=1 !initialization for lower bound
            do i_dist=1,N_dist
                !loop over distances (nearest, next nearest,... neighbor)
                Hmag=-io%trip(i_attrip)%val(i_dist)  !flip sign corresponding to previous implementation
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

                        if (present(neighbor_pos_list)) neighbor_pos_list(:,i_trip)=neigh%diff_vec(:,i_trip)

                        Htmp(offset_mag(1)+1,offset_mag(2)+2)= DMI(3) * Hmag
                        Htmp(offset_mag(1)+3,offset_mag(2)+1)= DMI(2) * Hmag
                        Htmp(offset_mag(1)+2,offset_mag(2)+3)= DMI(1) * Hmag

                        Htmp(offset_mag(1)+2,offset_mag(2)+1)=-DMI(3) * Hmag
                        Htmp(offset_mag(1)+1,offset_mag(2)+3)=-DMI(2) * Hmag
                        Htmp(offset_mag(1)+3,offset_mag(2)+2)=-DMI(1) * Hmag

                        Call get_coo(Htmp,val_tmp,ind_tmp)

                        !fill Hamiltonian type
                        Call Ham_tmp%init_connect(neigh%pairs(:,connect_bnd(1):connect_bnd(2)),val_tmp,ind_tmp,"MM",lat,2)
                        deallocate(val_tmp,ind_tmp)
                        Call Ham%add(Ham_tmp)

                        if (present(Ham_shell_pos)) then
                           Ham_shell_pos(:,:,i_trip)=Htmp
                        else
                           Call Ham_tmp%destroy()
                        endif
                    endif
                    connect_bnd(1)=connect_bnd(2)+1
                enddo 
            enddo
        enddo

        Ham%desc=ham_desc
        if(.not.Ham%is_set())then
            write(*,'(/A)') "DID NOT FIND A SINGLE DMI CONTRIBUTION, LEAVE THE DMI OUT IF YOU EXPECT THE DMI TO VANISH OR FIND A MISTAKE IN THIS CALCULATION"
            !continuing will give problems since the t_H-type Ham is not set but most probably will be assumed to be set in routine calling this subroutine
            STOP "CHECK INPUT"
        endif
        !set modes
        Call mode_set_rank1(Ham%mode_l,lat,"M")
        Call mode_set_rank1(Ham%mode_r,lat,"M")
    endif
end subroutine 


subroutine get_exchange_D_fft(H_fft,io,lat)
    !get heisenberg anti-symmetric exchange in fft_H Hamiltonian format
    !Since the anisotropy is localized localized in normal space this makes absolutely no sense unless used in combination with a delocalized Hamiltonian (dipolar-interaction) so that the evaluation is for free 
    use m_fft_H_base, only: fft_H
    use m_derived_types, only: lattice
    use m_setH_util, only: get_coo
    use m_neighbor_type, only: neighbors

    class(fft_H),intent(inout)  :: H_fft 
    type(io_H_D),intent(in)     :: io
    type(lattice),intent(in)    :: lat

    !fft parameters
    integer         :: Nmag             !number of magnetic atoms per unit-cell
    logical         :: period(3)        !consider as periodic or open boundary condition along each direction (T:period, F:open)
                                        ! (dim_lat(i)=1->period(i)=T, since the calculation in the periodic case is easier, but choice of supercell_vec still does not consider periodicity)
    integer         :: N_rep(3)         !number of states in each direction in the fourier transformation
    integer         :: Nk_tot           !number of state considered in FT (product of N_rep)
    integer         :: Kbd(2,3)         !boundaries of the K-operator

    real(8),allocatable :: Karr(:,:,:)  !K-operator to be FT'd (1:3*Nmag,1:3*Nmag,1:Nk_tot)


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

    integer         :: ind              !index in Karr space for given cell difference
    integer         :: ind3(3)          !index to calculate Karr position in each dimension
    integer         :: ind_mult(3)      !constant multiplicator to calculate ind from ind3 
    integer         :: i


    if(io%is_set)then
        !set some initial parameters locally for convencience
        write(output_unit,'(/2A)') "Start setting fft-Hamiltonian: ", ham_desc
        Nmag=lat%nmag
        period=lat%periodic.or.lat%dim_lat==1

        !set shape-dependent quantities of fft_H and get Kdb,N_rep
        Call H_fft%init_shape(3*lat%nmag,period,lat%dim_lat,Kbd,N_rep)
        Nk_tot=product(N_rep)

        !set local Hamiltonian 
        allocate(Karr(3*Nmag,3*Nmag,Nk_tot),source=0.0d0)
        ind_mult=[(product(N_rep(:i-1)),i=1,3)]
        N_attrip=size(io%trip)

        do i_attrip=1,N_attrip
            !loop over different connected atom types
            Call neigh%get(io%trip(i_attrip)%attype(1:2),io%trip(i_attrip)%dist,lat)
            !write information out
            Call io%trip(i_attrip)%prt(output_unit,'2X')
            Call neigh%prt(output_unit,'2X')
            !write information out
            Call io%trip(i_attrip)%prt(output_unit,'2X')
            Call neigh%prt(output_unit,'2X')
            N_dist=size(io%trip(i_attrip)%dist)
            N_dist=size(io%trip(i_attrip)%dist)
            i_trip=0
            connect_bnd=1 !initialization for lower bound
            do i_dist=1,N_dist
                !loop over distances (nearest, next nearest,... neighbor)
                Hmag=-io%trip(i_attrip)%val(i_dist)  !flip sign corresponding to previous implementation
                do i_shell=1,neigh%Nshell(i_dist)
                    !loop over all different connections with the same distance
                    i_trip=i_trip+1
                    Call get_DMI(neigh%at_pair(:,i_trip),neigh%pairs(:,connect_bnd(1)),io%trip(i_attrip)%attype(3),lat,DMI)
                    Call print_DMI(DMI,neigh%at_pair(:,i_trip),neigh%pairs(:,connect_bnd(1)),io%trip(i_attrip)%attype(3),io%trip(i_attrip)%dist(i_dist),lat)
                    connect_bnd(2)=neigh%ishell(i_trip)
                    if(norm2(DMI)>1.0d-8)then 
                                !    !set local Hamiltonian in basis of magnetic orderparameter
                        !find out which index in the Karr this entry corresponds to
                        atind_mag(1)=lat%cell%ind_mag(neigh%at_pair(1,i_trip))
                        atind_mag(2)=lat%cell%ind_mag(neigh%at_pair(2,i_trip))
                        offset_mag=(atind_mag-1)*3    !offset for magnetic index 
                        ind3=neigh%diff_cell(:,i_trip)
                        ind3=ind3-N_rep*floor(real(ind3,8)/lat%dim_lat)
                        ind3=ind3*ind_mult
                        ind=1+sum(ind3)

                        Karr(offset_mag(1)+1,offset_mag(2)+2,ind)= DMI(3) * Hmag
                        Karr(offset_mag(1)+3,offset_mag(2)+1,ind)= DMI(2) * Hmag
                        Karr(offset_mag(1)+2,offset_mag(2)+3,ind)= DMI(1) * Hmag

                        Karr(offset_mag(1)+2,offset_mag(2)+1,ind)=-DMI(3) * Hmag
                        Karr(offset_mag(1)+1,offset_mag(2)+3,ind)=-DMI(2) * Hmag
                        Karr(offset_mag(1)+3,offset_mag(2)+2,ind)=-DMI(1) * Hmag
                    endif
                    connect_bnd(1)=connect_bnd(2)+1
                enddo 
            enddo
        enddo
        Call H_fft%init_op(3*Nmag,Karr,ham_desc)
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
    Call get_neigh_distances(reshape(pos_center_uc,[3,1]),atpos_nonM,[1],lat,pairs,Nshell,distance)

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
