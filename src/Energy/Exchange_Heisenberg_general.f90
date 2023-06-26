module m_exchange_heisenberg_general
use m_input_H_types, only : io_H_Exchten
use m_rotation_matrix, only : rotate_matrix,check_rotate_matrix,rotation_matrix_real
use m_rotation, only : rotation_axis
use m_vector, only : norm
use m_sym_public
use m_symmetry_base
use m_sym_utils
use, intrinsic :: iso_fortran_env, only : output_unit
implicit none
private
public read_ExchG_input, get_exchange_ExchG, get_exchange_ExchG_fft

character(len=*),parameter  :: ham_desc="general magnetic exchange matrix"

contains

subroutine read_ExchG_input(io_param,fname,io)
    use m_io_utils
    integer,intent(in)              :: io_param
    character(len=*), intent(in)    :: fname
    type(io_H_Exchten),intent(out)  :: io

    Call get_parameter(io_param,fname,'magnetic_r2_tensor',io%pair,io%is_set)
    Call get_parameter(io_param,fname,'magnetic_tensor_fft',io%fft)
    Call get_parameter(io_param,fname,'c_H_Exchten',io%c_H_Exchten)

end subroutine

subroutine get_exchange_ExchG(Ham,io,lat,Ham_shell_pos,neighbor_pos_list)
    !get heisenberg symmetric exchange in t_H Hamiltonian format
    use m_H_public
    use m_derived_types, only: lattice
    use m_mode_public
    use m_setH_util, only: get_coo
    use m_neighbor_type, only: neighbors
    use m_constants, only : pi

    class(t_H),intent(inout)                       :: Ham  !Hamiltonian in which all contributions are added up
    type(io_H_Exchten),intent(in)                  :: io
    type(lattice),intent(in)                       :: lat
    real(8),optional,allocatable,intent(inout)     :: neighbor_pos_list(:,:)
    real(8),optional,allocatable,intent(inout)     :: Ham_shell_pos(:,:,:)

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
    real(8)         :: J(3,3),J_init(3,3)           !magnitude of Hamiltonian parameter
    integer         :: atind_mag(2)     !index of considered atom in basis of magnetic atoms (1:Nmag)
    integer         :: offset_mag(2)    !offset for start in dim_mode of chosed magnetic atom
    real(8)         :: bound_input(3)         ! first vector along which the interactions are given. It MUST be x=(1.0,0.0,0.0)
    real(8)         :: axis(3),angle,vec_tmp(3),scalaire,symop(3,3),chirality,periodic(3)
    integer         :: shell,k,i_op,i
    logical         :: found_sym
    class(pt_grp),allocatable :: my_symmetries
    character(len=30) :: name_sym

    if(io%is_set)then
        write(output_unit,'(/2A)') "Start setting Hamiltonian: ", ham_desc
        Call get_Htype(Ham_tmp)
        N_atpair=size(io%pair)

        allocate(Htmp(lat%M%dim_mode,lat%M%dim_mode))!local Hamiltonian modified for each shell/neighbor

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
          allocate(Ham_shell_pos(lat%M%dim_mode,lat%M%dim_mode,i_pair))
          Ham_shell_pos=0.0d0
          allocate(neighbor_pos_list(3,i_pair))
        endif

        !convert the periodic boundary conditions into reals for the chirality check
        periodic=0.0d0
        do i=1,3
           if (.not.lat%periodic(i)) periodic(i)=1.0d0
        enddo

        ! get the symmetries
        call set_sym_type(my_symmetries)
        call my_symmetries%read_sym('symmetries.out')

        do i_atpair=1,N_atpair
            !loop over different connected atom types
            !get neighbors
            Call neigh%get(io%pair(i_atpair)%attype,io%pair(i_atpair)%dist,lat)

            !write information out
            Call io%pair(i_atpair)%prt(output_unit,'2X')
            Call neigh%prt(output_unit,'2X')
            write(output_unit,*)
            !initialize variables
            N_dist=size(io%pair(i_atpair)%dist)
            i_pair=0
            connect_bnd=1 !initialization for lower bound

            do i_dist=1,N_dist
                !loop over distances (nearest, next nearest,... neighbor
                J_init=reshape(io%pair(i_atpair)%val(:,i_dist),(/3,3/))*io%c_H_Exchten
                bound_input=io%pair(i_atpair)%bound(:,i_dist)/norm(io%pair(i_atpair)%bound(:,i_dist))
                shell=0

                do i_shell=1,neigh%Nshell(i_dist)
                    !loop over all different connections with the same distance
                    i_pair=i_pair+1
                    shell=1+shell

                    !set local Hamiltonian in basis of magnetic orderparameter
                    atind_mag(1)=lat%cell%ind_mag(neigh%at_pair(1,i_pair))
                    atind_mag(2)=lat%cell%ind_mag(neigh%at_pair(2,i_pair))
                    offset_mag=(atind_mag-1)*3    !offset for magnetic index
                    Htmp=0.0d0

                    if (present(neighbor_pos_list)) neighbor_pos_list(:,i_pair)=neigh%diff_vec(:,i_pair)

                    ! rotate the exchange matrix to align it with the neighbor direction

                    vec_tmp=neigh%diff_vec(:,i_pair)/norm(neigh%diff_vec(:,i_pair))

                    do k=1,my_symmetries%n_sym
                       call check_rotate_matrix(my_symmetries%rotmat(k)%mat,bound_input,vec_tmp,found_sym,my_symmetries%tol_sym)
                       if (found_sym) then
                          i_op=k
                          exit
                         else
                          if (k.eq.my_symmetries%n_sym) STOP 'symmetry operation not found in Exchange_Heisenberg_general'
                       endif
                    enddo

                    symop=my_symmetries%rotmat(i_op)%mat
                    name_sym=my_symmetries%rotmat(i_op)%name
                    chirality=dot_product(matmul(symop,periodic),periodic)

                    call rotate_exchange(J,J_init,symop,chirality)

                    write(output_unit,'(A,I6,A)')   ' Applying exchange tensor along bound ',i_shell,':'
                    write(output_unit,'(2A)')       ' Applying symmetry operation ', trim(name_sym)
                    write(output_unit,'(A,2I6)')    '  atom types:', neigh%at_pair(:,shell)
                    write(output_unit,'(A,2I6)')    '  distance  :', io%pair(i_atpair)%dist(i_dist)
                    write(output_unit,'(A,9E16.8)') '  energy    :', J
                    write(output_unit,'(A,3E16.8/)') ' along the bound :', neigh%diff_vec(:,i_pair)

                    Htmp(offset_mag(1)+1:offset_mag(1)+3,offset_mag(2)+1:offset_mag(2)+3)=J
                    connect_bnd(2)=neigh%ishell(i_pair)

                    Call get_coo(Htmp,val_tmp,ind_tmp)

                    !fill Hamiltonian type
                    Call Ham_tmp%init_connect(neigh%pairs(:,connect_bnd(1):connect_bnd(2)),val_tmp,ind_tmp,"MM",lat,2)
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
        Call mode_set_rank1(Ham%mode_l,lat,"M")
        Call mode_set_rank1(Ham%mode_r,lat,"M")
    endif

end subroutine

subroutine get_exchange_ExchG_fft(H_fft,io,lat)
    !get heisenberg symmetric exchange in fft_H Hamiltonian format
    !Since the anisotropy is localized localized in normal space this makes absolutely no sense unless used in combination with a delocalized Hamiltonian (dipolar-interaction) so that the evaluation is for free
    use m_fft_H_base, only: fft_H
    use m_derived_types, only: lattice
    use m_setH_util, only: get_coo
    use m_neighbor_type, only: neighbors
    use m_constants, only : pi

    class(fft_H),intent(inout)    :: H_fft
    type(io_H_Exchten),intent(in) :: io
    type(lattice),intent(in)      :: lat

    !fft parameters
    integer         :: Nmag             !number of magnetic atoms per unit-cell
    logical         :: period(3)        !consider as periodic or open boundary condition along each direction (T:period, F:open)
                                        ! (dim_lat(i)=1->period(i)=T, since the calculation in the periodic case is easier, but choice of supercell_vec still does not consider periodicity)
    integer         :: N_rep(3)         !number of states in each direction in the fourier transformation
    integer         :: Nk_tot           !number of state considered in FT (product of N_rep)
    integer         :: Kbd(2,3)         !boundaries of the K-operator

    real(8),allocatable :: Karr(:,:,:)  !K-operator to be FT'd (1:3*Nmag,1:3*Nmag,1:Nk_tot)

    integer         :: i_atpair,N_atpair    !loop parameters which atom-type connection are considered (different neighbor types)
    integer         :: i_dist,N_dist        !loop parameters which  connection are considered (different neighbor types)
    integer         :: i_pair           !loop keeping track which unique connection between the same atom types is considered (indexes "number shells" in neighbors-type)
    integer         :: i_shell          !counting the number of unique connection for given atom types and a distance
    integer         :: connect_bnd(2)   !indices keeping track of which pairs are used for the particular connection
    type(neighbors) :: neigh            !all neighbor information for a given atom-type pair
    real(8)         :: J(3,3),J_init(3,3)           !magnitude of Hamiltonian parameter
    integer         :: atind_mag(2)     !index of considered atom in basis of magnetic atoms (1:Nmag)
    integer         :: offset_mag(2)    !offset for start in dim_mode of chosed magnetic atom

    integer         :: ind              !index in Karr space for given cell difference
    integer         :: ind3(3)          !index to calculate Karr position in each dimension
    integer         :: ind_mult(3)      !constant multiplicator to calculate ind from ind3
    integer         :: i,k,i_op,shell
    real(8)         :: axis(3),angle,vec_tmp(3),bound_input(3),scalaire,symop(3,3),chirality,periodic(3)
    logical         :: found_sym
    class(pt_grp),allocatable :: my_symmetries
    character(len=30) :: name_sym

    if(io%is_set)then
        write(output_unit,'(/2A)') "Start setting fft-Hamiltonian: ", ham_desc
        !set some initial parameters locally for convencience
        Nmag=lat%nmag
        call H_fft%set_mode_id(1,"M")
        period=lat%periodic.or.lat%dim_lat==1

        !convert the periodic boundary conditions into reals for the chirality check
        periodic=0.0d0
        do i=1,3
           if (.not.lat%periodic(i)) periodic(i)=1.0d0
        enddo

        ! get the symmetries
        call set_sym_type(my_symmetries)
        call my_symmetries%read_sym('symmetries.out')

        !set shape-dependent quantities of fft_H and get Kdb,N_rep
        Call H_fft%init_shape("M",3*lat%nmag,period,lat%dim_lat,Kbd,N_rep)
        Nk_tot=product(N_rep)

        !set local Hamiltonian
        allocate(Karr(3*Nmag,3*Nmag,Nk_tot),source=0.0d0)
        ind_mult=[(product(N_rep(:i-1)),i=1,3)]
        N_atpair=size(io%pair)
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
                J_init=io%c_H_Exchten**reshape(io%pair(i_atpair)%val(:,i_dist),(/3,3/)) !negative sign for compatibility with previous implementatoins
                bound_input=io%pair(i_atpair)%bound(:,i_dist)/norm(io%pair(i_atpair)%bound(:,i_dist))
                shell=0

                do i_shell=1,neigh%Nshell(i_dist)
                    !loop over all different connections with the same distance
                    i_pair=i_pair+1
                    shell=1+shell

                    !find out which index in the Karr this entry corresponds to
                    atind_mag(1)=lat%cell%ind_mag(neigh%at_pair(1,i_pair))
                    atind_mag(2)=lat%cell%ind_mag(neigh%at_pair(2,i_pair))
                    offset_mag=(atind_mag-1)*3    !offset for magnetic index
                    ind3=neigh%diff_cell(:,i_pair)
                    ind3=ind3-N_rep*floor(real(ind3,8)/lat%dim_lat)
                    ind3=ind3*ind_mult
                    ind=1+sum(ind3)

                    ! rotate the exchange matrix to align it with the neighbor direction
                     ! rotation axis
                    vec_tmp=neigh%diff_vec(:,i_pair)/norm(neigh%diff_vec(:,i_pair))

                    do k=1,my_symmetries%n_sym
                       call check_rotate_matrix(my_symmetries%rotmat(k)%mat,bound_input,vec_tmp,found_sym,my_symmetries%tol_sym)
                       if (found_sym) then
                          i_op=k
                          exit
                         else
                          if (k.eq.my_symmetries%n_sym) STOP 'symmetry operation not found in Exchange_Heisenberg_general'
                       endif
                    enddo

                    symop=my_symmetries%rotmat(i_op)%mat
                    name_sym=my_symmetries%rotmat(i_op)%name
                    chirality=dot_product(matmul(symop,periodic),periodic)

                    call rotate_exchange(J,J_init,symop,chirality)

                    write(output_unit,'(A,I6,A)')   ' Applying exchange tensor along bound ',i_shell,':'
                    write(output_unit,'(2A)')       ' Applying symmetry operation ', trim(name_sym)
                    write(output_unit,'(A,2I6)')    '  atom types:', neigh%at_pair(:,shell)
                    write(output_unit,'(A,2I6)')    '  distance  :', io%pair(i_atpair)%dist(i_dist)
                    write(output_unit,'(A,9E16.8)') '  energy    :', J
                    write(output_unit,'(A,3E16.8/)') ' along the bound :', neigh%diff_vec(:,i_pair)



                    !set the contributions in the operator array
                    Karr(offset_mag(1)+1:offset_mag(1)+3,offset_mag(2)+1:offset_mag(2)+3,ind)=J
                enddo
            enddo
        enddo
        Call H_fft%init_op(3*Nmag,Karr,ham_desc)
    endif
end subroutine

end module
