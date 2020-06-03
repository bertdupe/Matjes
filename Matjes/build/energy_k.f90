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
            complex(kind=16), allocatable :: W(:), VL(:,:), VR(:,:), WORK(:)

            N = size(all_E_k, 1)
            allocate( W(N), VL(N,N), VR(N,N), WORK(4*N) )

            ! Before diagonalising the Hamiltonian, we first have to Fourier transform it
            call CGEEV('N', 'V', N, get_FFT(all_E_k, pos, kvector_pos, dim_mode, sense), LDA, W, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO)

            deallocate( W, VL, VR, WORK )
        end subroutine diagonalise_H_k



        ! Function implementing the FD distribution for a given energy
        real(kind=8) function fermi_distrib(E, energy, kt)
            implicit none
            real(kind=8) :: E(:)
            real(kind=8) :: kt
            real(kind=8) :: energy, fermi_distrib

            fermi_distrib = 1.0/( 1.0 + exp((E-energy)/kt )
        end function fermi_distrib



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
        subroutine compute_Etot( Etot, E, eps_nk , kt)
            implicit none
            real(kind=8), intent(in) :: eps_nk(:), kt
            real(kind=8), intent(inout) :: Etot

            ! Internal variable
            integer :: i

            do i=1, size(eps_nk)
                Etot = Etot + eps_nk(i)*fermi_distrib(E, eps_nk(i), kt)
            enddo

        end subroutine compute_Etot

end module m_energy_k
