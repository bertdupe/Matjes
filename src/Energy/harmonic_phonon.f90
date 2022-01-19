module m_harmonic_phonon
use, intrinsic :: iso_fortran_env, only : output_unit
use m_input_H_types, only: io_H_Ph
implicit none

character(len=*),parameter  :: ham_desc="harmonic phonon"
logical :: read_from_file = .False.
character(len=30), parameter :: fname_phonon = 'phonon_harmonic.in'

private
public :: read_F_input, get_Forces_F

contains

subroutine read_F_input(io_param,fname,io)
    use m_io_utils
    integer,intent(in)              :: io_param
    character(len=*), intent(in)    :: fname
    type(io_H_Ph),intent(out)        :: io

    write(6,'(a)') 'reading forces in Ha/Bohr^2'
    Call get_parameter(io_param,fname,'phonon_harmonic',io%pair,io%is_set)
    if (io%is_set) Call get_parameter(io_param,fname,'c_ph',io%c_ph)

    inquire(file=fname_phonon,exist=read_from_file)
    if (read_from_file) write(6,'(a)') 'reading phonon from phonon_harmonic.in'

end subroutine

subroutine get_Forces_F(Ham,io,lat,Ham_shell_pos,neighbor_pos_list)
    !get coupling in t_H Hamiltonian format
    use m_H_public
    use m_derived_types
    use m_setH_util,only: get_coo
    use m_neighbor_type, only: neighbors
    use m_forces_from_file, only: get_forces_file
    use m_mode_public
    use m_vector, only : norm

    class(t_H),intent(inout)                         :: Ham
    type(io_H_Ph),intent(in)                         :: io
    type(lattice),intent(in)                         :: lat
    real(8),optional,allocatable,intent(inout)       :: neighbor_pos_list(:,:)
    real(8),optional,allocatable,intent(inout)       :: Ham_shell_pos(:,:,:)

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
    real(8)         :: F                !magnitude of Hamiltonian parameter
    integer         :: atind_ph(2)      !index of considered atom in basis of phonon atoms (1:NPh)
    integer         :: offset_ph(2)     !offset for start in dim_mode of chosed phonon atom
    real(8)         :: norm_vec_neigh,vec_neigh(3)

    ! conversion factor Ha/Bohr2 to eV/nm2
    ! 1 Ha/Bohr = 51.42208619083232 eV/Angstrom
    real(8), parameter  :: HaBohrsq_to_Evnmsq = 9717.38d0
    real(8), parameter  :: HaBohr_to_Evnm = 514.2208619083232d0

    if(io%is_set)then
        write(output_unit,'(/2A)') "Start setting Hamiltonian: ", ham_desc
        Call get_Htype(Ham_tmp)
        N_atpair=size(io%pair)
        allocate(Htmp(lat%u%dim_mode,lat%u%dim_mode))!local Hamiltonian modified for each shell/neighbor

        if (present(Ham_shell_pos)) then
          write(output_unit,'(/2A)') "Preparing the Fourier Transform of Hamiltonian: ", ham_desc
          i_pair=0
          do i_atpair=1,N_atpair
             Call neigh%get(io%pair(i_atpair)%attype,io%pair(i_atpair)%dist,lat)
             N_dist=size(io%pair(i_atpair)%dist)
             do i_dist=1,N_dist
                do i_shell=1,neigh%Nshell(i_dist)
                   i_pair=i_pair+1
                enddo
             enddo
          enddo
          allocate(Ham_shell_pos(lat%u%dim_mode,lat%u%dim_mode,i_pair))
          Ham_shell_pos=0.0d0
          allocate(neighbor_pos_list(3,i_pair))
        endif

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
                F=io%pair(i_atpair)%val(i_dist)

                do i_shell=1,neigh%Nshell(i_dist)
                    !loop over all different connections with the same distance
                    i_pair=i_pair+1

                    !set local Hamiltonian in basis of displacement orderparameter
                    atind_ph(1)=lat%cell%ind_ph(neigh%at_pair(1,i_pair))
                    atind_ph(2)=lat%cell%ind_ph(neigh%at_pair(2,i_pair))
                    vec_neigh=neigh%diff_vec(:,i_pair)

                    Htmp=0.0d0
                    offset_ph=(atind_ph-1)*3
                    if (present(neighbor_pos_list)) neighbor_pos_list(:,i_pair)=vec_neigh
                    if (read_from_file) then
                       ! (Hamiltonian , name of the file , relative position of the neighbor , offset)
                       write(output_unit,'(/2A)') 'Hamiltonian from DFT is in Ha/Bohr2'
                       F=F*HaBohrsq_to_Evnmsq
                       call get_forces_file(Htmp,fname_phonon,neigh%diff_vec(:,i_pair),offset_ph)
                       Htmp=Htmp*HaBohrsq_to_Evnmsq
                    else
                       write(output_unit,'(/2A)') 'Hamiltonian directly from input is in eV/nm2'
                       norm_vec_neigh=norm(vec_neigh)
                       Htmp(offset_ph(1)+1,offset_ph(2)+1)=F*abs(vec_neigh(1))/norm_vec_neigh
                       Htmp(offset_ph(1)+2,offset_ph(2)+2)=F*abs(vec_neigh(2))/norm_vec_neigh
                       Htmp(offset_ph(1)+3,offset_ph(2)+3)=F*abs(vec_neigh(3))/norm_vec_neigh
                    endif
                    connect_bnd(2)=neigh%ishell(i_pair)

                    Htmp=Htmp*io%c_ph


                    Call get_coo(Htmp,val_tmp,ind_tmp)

                    !fill Hamiltonian type
                    Call Ham_tmp%init_connect(neigh%pairs(:,connect_bnd(1):connect_bnd(2)),val_tmp,ind_tmp,"UU",lat,0)
                    deallocate(val_tmp,ind_tmp)
                    Call Ham%add(Ham_tmp)
                    if (present(Ham_shell_pos)) then
                       Ham_shell_pos(:,:,i_pair)=Htmp
                    else
                       Call Ham_tmp%destroy()
                    endif
                    connect_bnd(1)=connect_bnd(2)+1
                enddo 
            enddo
        enddo
        Ham%desc=ham_desc
        !set modes
        Call mode_set_rank1(Ham%mode_l,lat,"U")
        Call mode_set_rank1(Ham%mode_r,lat,"U")
    endif

end subroutine

end module m_harmonic_phonon
