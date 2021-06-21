module m_ASR_phonon
use m_input_H_types, only: io_U_ASR
use m_forces_from_file, only: get_ASR_file

implicit none
logical :: read_from_file = .False.

private :: read_from_file
public :: get_ASR_Ph,read_ASR_Ph_input

contains

subroutine read_ASR_Ph_input(io_param,fname,io)
    use m_io_utils
    integer,intent(in)              :: io_param
    character(len=*), intent(in)    :: fname
    type(io_U_ASR),intent(out)      :: io
    logical                         :: cancel_ASR

    cancel_ASR=.False.
    Call get_parameter(io_param,fname,'cancel_ASR',cancel_ASR)
    if (.not.cancel_ASR) then
      Call get_parameter(io_param,fname,'phonon_harmonic',io%pair,io%is_set)

      inquire(file='phonon_harmonic.in',exist=read_from_file)
      if (read_from_file) write(6,'(a)') 'reading ASR from phonon_harmonic.in'

    endif

end subroutine

subroutine get_ASR_Ph(Ham,io,lat)
    !get anisotropy in t_H Hamiltonian format
    use m_H_public
    use m_derived_types, only: lattice
    use m_setH_util, only: get_coo
    use m_neighbor_type, only: neighbors
    use m_mode_public

    class(t_H),intent(inout)    :: Ham
    type(io_U_ASR),intent(in)   :: io
    type(lattice),intent(in)    :: lat
    !local

    class(t_H),allocatable    :: Ham_tmp    !temporary Hamiltonian type used to add up Ham

    integer :: i_atpair,i_dist,N_dist,i_shell,i_pair
    real(8),allocatable :: Htmp(:,:)
    real(8),allocatable :: val_tmp(:)
    integer,allocatable :: ind_tmp(:,:)
    integer             :: N_atpair    ! nb of shells that have to be taken into ac count
    type(neighbors)     :: neigh            !all neighbor information for a given atom-type pair
    real(8)             :: F                !magnitude of Hamiltonian parameter
    integer             :: connect_bnd(2)   !indices keeping track of which pairs are used for the particular connection
    integer             :: atind_ph(2)     !index of considered atom in basis of magnetic atoms (1:Nmag)
    integer,allocatable :: all_pairs(:,:)
    integer             :: offset_ph(2)     !offset for start in dim_mode of chosed phonon atom

    ! conversion factor Ha/Bohr2 to eV/nm2
    ! 1 Ha/Bohr = 51.42208619083232 eV/Angstrom
    real(8), parameter  :: HaBohrsq_to_Evnmsq = 9717.38d0
    real(8), parameter  :: HaBohr_to_Evnm = 514.2208619083232d0


    if(io%is_set)then
        Call get_Htype(Ham_tmp)
        N_atpair=size(io%pair)
        !set local Hamiltonian
        allocate(Htmp(lat%u%dim_mode,lat%u%dim_mode),source=0.d0)

        if (read_from_file) then

            do i_atpair=1,N_atpair
                !loop over different connected atom types
                Call neigh%get(io%pair(i_atpair)%attype,io%pair(i_atpair)%dist,lat)
                N_dist=size(io%pair(i_atpair)%dist)
                connect_bnd=1 !initialization for lower bound
                connect_bnd(2)=neigh%ishell(i_atpair)

                !set local Hamiltonian in basis of magnetic orderparameter
                atind_ph(1)=lat%cell%ind_ph(neigh%at_pair(1,i_atpair))
                atind_ph(2)=lat%cell%ind_ph(neigh%at_pair(2,i_atpair))
                offset_ph=(atind_ph-1)*3

                allocate( all_pairs,source=neigh%pairs(:,connect_bnd(1):connect_bnd(2)) )
                Htmp=0.0d0
                all_pairs(2,:)=neigh%pairs(1,connect_bnd(1):connect_bnd(2))

                call get_ASR_file(Htmp,'phonon_harmonic.in',offset_ph)
                Htmp=Htmp*HaBohrsq_to_Evnmsq
                Call get_coo(Htmp,val_tmp,ind_tmp)

                !fill Hamiltonian type
                Call Ham_tmp%init_connect(all_pairs,val_tmp,ind_tmp,"UU",lat,2)
                deallocate(val_tmp,ind_tmp)
                Call Ham%add(Ham_tmp)
                Call Ham_tmp%destroy()

                deallocate( all_pairs )
            enddo

        else
          do i_atpair=1,N_atpair
            !loop over different connected atom types
            Call neigh%get(io%pair(i_atpair)%attype,io%pair(i_atpair)%dist,lat)
            N_dist=size(io%pair(i_atpair)%dist)
            connect_bnd=1 !initialization for lower bound
            i_pair=0

            do i_dist=1,N_dist
                !loop over distances (nearest, next nearest,... neighbor) also called shell
                F=io%pair(i_atpair)%val(i_dist)*HaBohrsq_to_Evnmsq

                do i_shell=1,neigh%Nshell(i_dist)
                    write(*,*) 'i_shell',i_shell
                    !loop over all different connections with the same distance
                    i_pair=i_pair+1
                    connect_bnd(2)=neigh%ishell(i_pair)

                    !set local Hamiltonian in basis of magnetic orderparameter
                    atind_ph(1)=lat%cell%ind_ph(neigh%at_pair(1,i_pair))
                    atind_ph(2)=lat%cell%ind_ph(neigh%at_pair(2,i_pair))

                    allocate( all_pairs,source=neigh%pairs(:,connect_bnd(1):connect_bnd(2)) )
                    Htmp=0.0d0

                    Htmp(atind_ph(1)*3-2,atind_ph(1)*3-2)=io%c_ASR*F/2.0d0
                    Htmp(atind_ph(1)*3-1,atind_ph(1)*3-1)=io%c_ASR*F/2.0d0
                    Htmp(atind_ph(1)*3  ,atind_ph(1)*3  )=io%c_ASR*F/2.0d0

                    all_pairs(2,:)=neigh%pairs(1,connect_bnd(1):connect_bnd(2))

                    Call get_coo(Htmp,val_tmp,ind_tmp)

                    !fill Hamiltonian type
                    Call Ham_tmp%init_connect(all_pairs,val_tmp,ind_tmp,"UU",lat,2)
                    deallocate(val_tmp,ind_tmp)
                    Call Ham%add(Ham_tmp)
                    Call Ham_tmp%destroy()


                    all_pairs(1,:)=neigh%pairs(2,connect_bnd(1):connect_bnd(2))
                    all_pairs(2,:)=neigh%pairs(2,connect_bnd(1):connect_bnd(2))

                    Htmp=0.0d0

                    Htmp(atind_ph(2)*3-2,atind_ph(2)*3-2)=io%c_ASR*F/2.0d0
                    Htmp(atind_ph(2)*3-1,atind_ph(2)*3-1)=io%c_ASR*F/2.0d0
                    Htmp(atind_ph(2)*3  ,atind_ph(2)*3  )=io%c_ASR*F/2.0d0

                    Call get_coo(Htmp,val_tmp,ind_tmp)

                    !fill Hamiltonian type
                    Call Ham_tmp%init_connect(all_pairs,val_tmp,ind_tmp,"UU",lat,2)
                    deallocate(val_tmp,ind_tmp)
                    Call Ham%add(Ham_tmp)
                    Call Ham_tmp%destroy()

                    connect_bnd(1)=connect_bnd(2)+1
                    deallocate( all_pairs )

                enddo
            enddo
          enddo
        endif
        Ham%desc="ASR phonon"
        !set modes
        Call mode_set_rank1(Ham%mode_l,lat,"U")
        Call mode_set_rank1(Ham%mode_r,lat,"U")
    endif

end subroutine

end module m_ASR_phonon
