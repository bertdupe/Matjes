module m_energy_k
    use m_energy_commons, only : energy
    use m_fftw

    integer, allocatable, dimension(:,:) :: n_lines

    ! This matrix will contain all the Hamiltonians at the right places
    complex(kind=16), allocatable, dimension(:,:) :: all_E_k

    private
    public :: rewrite_H_k, diagonalise_H_k, fermi_distrib, compute_Etot
    contains
        subroutine rewrite_H_k(dim_mode)
            implicit none
            integer, intent(in) :: dim_mode

            ! Internal variables
            integer :: i, j, k, N_neighbours, N_cells

            N_neighbours = size( energy%line, 1 )
            N_cells = size( energy%line, 2 )
            allocate(all_E_k(dim_mode*N_cells, dim_mode*N_cells), n_lines(N_neighbours, N_cells))
            all_E_k = cmplx(0.0d0, kind=16)
            n_lines = energy%line !contains the columns of the non-zero Hamiltonians

            do i=1, N_cells
                do j=1, N_neighbours
                    k = energy%line(j, i) ! "k" is the index in the vector of the "j" non-zero Hamiltonian for site "i"
                    ! energy%value(1,k)%order_op(1)%Op_loc corresponds to the onsite
                    ! energy%value(2,k)%order_op(1)%Op_loc corresponds to the first neighbour
                    ! energy%value(3,k)%order_op(1)%Op_loc corresponds to the second neighbour
                    ! etc.
                    ! This means that the first argument in "value" corresponds to the shell
                    ! The matrix all_E_k will have the structure
                    !   [H0 H1 H2 0 .......... 0 H2 H1]
                    !   [H1 H0 H1 H2 0 .......... 0 H2]
                    !   [.............................]
                    !   [.............................]
                    !   [.............................]
                    !   [H1 H2 0 ...........0 H2 H1 H0]
!                    all_E_k( (k-1)*dim_mode+1:k*dim_mode, (i-1)*dim_mode+1:i*dim_mode ) = energy%value(j,k)%order_op(1)%Op_loc ! energy%value(j,k) is the "k" Hamiltonian (non-zero) in the "j" shell
                    all_E_k( (i-1)*dim_mode+1:i*dim_mode, (k-1)*dim_mode+1:k*dim_mode ) = energy%value(j,k)%order_op(1)%Op_loc ! energy%value(j,k) is the "k" Hamiltonian (non-zero) in the "j" shell
!write(*,*) 'all_E_k(', (i-1)*dim_mode+1, ':', i*dim_mode, ',', (k-1)*dim_mode+1, ':', k*dim_mode, ') = ', all_E_k( (i-1)*dim_mode+1:i*dim_mode, (k-1)*dim_mode+1:k*dim_mode )
!write(*,*) ''
!write(*,*) ''
!write(*,*) ''
                enddo
            enddo
        end subroutine rewrite_H_k



        subroutine diagonalise_H_k(kvector_pos, pos, dim_mode, sense)
            implicit none
            integer :: kvector_pos, dim_mode
            real(kind=8) :: sense
            real(kind=8), intent(in) :: pos(:,:)

            ! Internal variable
            integer :: i

            ! JOBVL: left eigenvec of A are computed ('V') or not ('N')
            ! JOBVR: right eigenvec of A are computed ('V') or not ('N')
            ! N: order of matrix A
            ! A: N-by-N input matrix. Is overwritten at output
            ! LDA: leading dimension of A
            ! W: eigenvalues of A
            ! VL: left eigenvec if JOVL is 'V'. Not referenced otherwise
            ! LDVL: leading dimension of array VL
            ! VR: right eigenvec if JOVL is 'V'. Not referenced otherwise
            ! LDVR: leading dimension of array VR
            ! WORK: ???
            ! LWORK: dimension of array WORK
            ! RWORK: real array of dimension 2*N
            ! INFO: computation information
            ! CGEEV(JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO)
            integer :: N, LDA, LDVL, LDVR, LWORK, RWORK, INFO
            complex(kind=16), allocatable :: W(:), VL(:,:), VR(:,:), WORK(:), H_complex(:,:)

            N = size(all_E_k, 1)
            allocate( W(N), VL(N,N), VR(N,N), WORK(4*N), H_complex(N,N))
            ! Before diagonalising the Hamiltonian, we first have to Fourier transform it
            H_complex=get_FFT(all_E_k, pos, kvector_pos, dim_mode, sense)

#ifdef CPP_BLAS
            ! diagonalising the Hamiltonian
            call CGEEV('N', 'V', N, H_complex , LDA, W, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO)
#endif

            deallocate( W, VL, VR, WORK )
        end subroutine diagonalise_H_k



        ! Subroutine computing the total energy contained in the system
        ! The total energy contained in the system is given by
        !   sum_{n,k} f_{FD}(epsilon_{n,k}) epsilon_{n,k}
        ! where epsilon_{n,k} are the eigenenergies and f_{FD} is the
        ! Fermi-Dirac distribution corresponding to that energy.
        ! Input:
        !   _ E is the energy vector
        !   _ eps_nk is the vector containing all the eigenvalues
        !   _ Etot is the total energy contained in the system
        ! Output:
        !   _ Etot is the total energy contained in the system
        subroutine compute_Etot( Etot, E, eps_nk , kt, N_electrons)
            implicit none
            complex(kind=16), intent(in) :: E(:)
            integer, intent(in) :: N_electrons
            real(kind=8), intent(in) :: eps_nk(:), kt
            real(kind=8), intent(out) :: Etot

            ! Internal variable
            integer :: i
            real(kind=8) :: fermi_level

!            call compute_Fermi_level(eps_nk, N_electrons, fermi_level, kt)

            Etot=0.0d0
            do i=1, size(eps_nk)
                Etot = Etot + eps_nk(i)*fermi_distrib(fermi_level, eps_nk(i), kt)
            enddo

        end subroutine compute_Etot

! Function computing the Fermi energy of the system
!        subroutine compute_Fermi_level(eps_nk, N_electrons, fermi_level, kt)
!            use m_sort
!            implicit none
!            real(kind=8), intent(out) :: fermi_level
!            complex(kind=16), intent(in) :: eps_nk(:)
!            integer, intent(in) :: N_electrons
!            ! Internal variable
!            integer :: i,size_ham
!            real(kind=8) :: precision_Ef, tmp_sum
!            real(kind=16), allocatable :: eps_nk_real(:)
!
!            size_ham=size(eps_nk)
!            allocate(eps_nk_real(size_ham))
!            eps_nk_real=real(eps_nk)
!            eps_nk_real=sort(eps_nk_real)
!
!            fermi_level = real(eps_nk_real(1),kind=8)
!            precision_Ef = 1.0d-3
!
!            tmp_sum = 0.0d0
!            while (abs(tmp_sum-real(N_electrons).)
!            do i=1, size_ham
!                tmp_sum = tmp_sum + fermi_distrib(fermi_level, eps_nk(i), kt)
!            enddo
!        end subroutine compute_Fermi_level






        ! Function implementing the FD distribution for a given energy
        real(kind=8) function fermi_distrib(E_F, energy, kt)
            implicit none
            real(kind=8) :: E_F
            real(kind=8) :: kt
            real(kind=8) :: energy

            fermi_distrib=0.0d0
            if (kt.lt.1.0d-8) then
              if (energy.lt.E_F) fermi_distrib=1.0d0
            else
              fermi_distrib = 1.0d0/( 1.0d0 + exp((E_F-energy)/kt ) )
            endif
        end function fermi_distrib

end module m_energy_k