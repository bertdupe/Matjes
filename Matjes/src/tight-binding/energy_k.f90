module m_energy_k
    use m_basic_types, only : vec_point
    use m_energy_commons, only : energy
    use m_fftw, only: get_FFT
    use m_J_sd_exchange
    use m_energy_set_real, only: set_Hr, get_Hr

    implicit none

    private
    public :: diagonalise_H_k_list,set_dist_neigh,get_energy_kpts !, fermi_distrib, compute_Etot
    contains

        subroutine get_energy_kpts(klist,dimH,tb_ext,pos,mode_mag,eigval)
            real(8),intent(in)                      ::  klist(:,:)
            integer,intent(in)                      ::  dimH
            integer,intent(in)                      ::  tb_ext(2)
            real(8),intent(in)                      ::  pos(:,:)
            type(vec_point),intent(in)              ::  mode_mag(dimH)
            real(8),intent(inout),allocatable       ::  eigval(:,:)

            complex(8)                  ::  Hr(dimH,dimH)
            integer                     ::  dim_mode
            
            !set real space Hamiltonian
            !if used more often, set this in advance
            Call set_Hr(dimH,tb_ext,mode_mag)
            Call get_Hr(dimH,Hr)
            !get actual eigenvalues
            allocate( eigval(dimH,size(klist,2)),source=0.0d0)
            dim_mode=TB_ext(2)-TB_ext(1)+1
            Call diagonalise_H_k_list(Hr,pos,klist,dim_mode,-1.0d0,eigval)
        end subroutine

        subroutine set_dist_neigh(dist_neigh,pos)
            real(8),intent(out),allocatable   ::  dist_neigh(:,:)
            real(8),intent(in)                  ::  pos(:,:)
            integer :: N_neighbours,i,j

            N_neighbours = size( energy%line, 1 )
            allocate(dist_neigh(3,N_neighbours),source=0.0d0)
            do i=1, N_neighbours
               j = energy%line(i, 1)
               dist_neigh(:,i) = pos(:,j)
            enddo
        end subroutine

#ifdef CPP_LAPACK
        subroutine diagonalise_H_k_list(Hr,pos,kpts,dim_mode,sense,eigval)
            external :: ZHEEV
            complex(8),intent(in)   ::  Hr(:,:)
            real(8),intent(in)      ::  pos(:,:)
            real(8),intent(in)      ::  kpts(:,:)
            integer, intent(in)     ::  dim_mode
            real(kind=8), intent(in) :: sense
            real(kind=8), intent(inout),allocatable :: eigval(:,:)

            integer :: i
            integer :: N,Nkpt, INFO
            complex(kind=8), allocatable :: WORK(:), H_complex(:,:)
            real(kind=8), allocatable :: RWORK(:)

            Nkpt=size(kpts,2)
            N = size(Hr, 1)
            if(.not. allocated(eigval)) allocate(eigval(N,Nkpt),source=0.0d0)
            if(size(eigval,2)/= Nkpt) STOP "eigval array does not have correct dimension(Nkpts)"
            allocate( H_complex(N,N), RWORK(max(1,3*N-2)), WORK(2*N))
            do i=1,Nkpt
                H_complex=get_FFT(Hr, pos, kpts(:,i), dim_mode, sense)
                call ZHEEV( 'N', 'U', N, H_complex, N, eigval(:,i), WORK, size(Work), RWORK, INFO )
            enddo
        end subroutine
#endif

!#ifdef CPP_INTERNAL
!I HAVE NO CLUE IF THIS WORKS, SO I COMMENT IT OUT RESTRUCTURING THIS MODULE
!        subroutine diagonalise_H_k(kvector_pos, dim_mode, sense, eigval)
!            use m_eigen_val_vec
!            use m_invert
!            implicit none
!            integer :: kvector_pos, dim_mode
!            real(kind=8) :: sense
!            complex(kind=8), intent(out) :: eigval(:)
!
!            ! Internal variable
!            integer :: i,N
!
!            complex(kind=8), allocatable ::  H_complex(:,:), U(:,:)
!
!            N = size(all_E_k, 1)
!            allocate(H_complex(N,N),U(N,N))
!            U=0.0d0
!
!            ! Before diagonalising the Hamiltonian, we first have to Fourier transform it
!            H_complex=get_FFT(all_E_k, all_positions, kvector_pos, dim_mode, sense)
!            
!            call Jacobi(1.0d-7,N,H_complex,N,eigval,U,N,0)
!
!        end subroutine diagonalise_H_k
!#endif
!
!
!        ! Subroutine computing the total energy contained in the system
!        ! The total energy contained in the system is given by
!        !   sum_{n,k} f_{FD}(epsilon_{n,k}) epsilon_{n,k}
!        ! where epsilon_{n,k} are the eigenenergies and f_{FD} is the
!        ! Fermi-Dirac distribution corresponding to that energy.
!        ! Input:
!        !   _ E is the energy vector
!        !   _ eps_nk is the vector containing all the eigenvalues
!        !   _ Etot is the total energy contained in the system
!        ! Output:
!        !   _ Etot is the total energy contained in the system
!        subroutine compute_Etot( Etot, E, eps_nk , kt, N_electrons)
!            implicit none
!            complex(kind=8), intent(in) :: E(:)
!            integer, intent(in) :: N_electrons
!            real(kind=8), intent(in) :: eps_nk(:), kt
!            real(kind=8), intent(out) :: Etot
!
!            ! Internal variable
!            integer :: i
!            real(kind=8) :: fermi_level
!
!!            call compute_Fermi_level(eps_nk, N_electrons, fermi_level, kt)
!
!            Etot=0.0d0
!            do i=1, size(eps_nk)
!                Etot = Etot + eps_nk(i)*fermi_distrib(fermi_level, eps_nk(i), kt)
!            enddo
!
!        end subroutine compute_Etot
!
!
!
!        ! Function implementing the FD distribution for a given energy
!        real(kind=8) function fermi_distrib(E_F, energy, kt)
!            implicit none
!            real(kind=8) :: E_F
!            real(kind=8) :: kt
!            real(kind=8) :: energy
!
!            fermi_distrib=0.0d0
!            if (kt.lt.1.0d-8) then
!                if (energy.le.E_F) fermi_distrib=1.0d0
!            else
!                fermi_distrib = 1.0d0/( 1.0d0 + exp((energy-E_F)/kt ) )
!            endif
!        end function fermi_distrib

end module m_energy_k
