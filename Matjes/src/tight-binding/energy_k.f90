module m_energy_k
    use m_energy_commons, only : energy
    use m_fftw

    integer, allocatable, dimension(:,:) :: n_lines
    
    ! This matrix will contain all the Hamiltonians at the right places
    complex(kind=8), allocatable, dimension(:,:) :: all_E_k
    real(kind=8), allocatable :: all_positions(:,:) ! array containing the coordinates of all r-r'

    private
    public :: rewrite_H_k, diagonalise_H_k, fermi_distrib, compute_Etot, compute_Fermi_level
    contains
        subroutine rewrite_H_k(dim_mode,TB_pos_start,TB_pos_end,pos)
            implicit none
            integer, intent(in) :: dim_mode,TB_pos_start,TB_pos_end
            real(kind=8), intent(in) :: pos(:,:)

            ! Internal variables
            integer :: i, j, k, N_neighbours, N_cells

            N_neighbours = size( energy%line, 1 )
            N_cells = size( energy%line, 2 )
            allocate(all_E_k(dim_mode*N_cells, dim_mode*N_cells),all_positions(3,N_neighbours))
            all_E_k = cmplx(0.0d0, kind=8)

            do j=1, N_neighbours
               k = energy%line(j, 1)
               all_positions(:,j) = pos(:,k)
            enddo

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
                    all_E_k( (i-1)*dim_mode+1:i*dim_mode, (k-1)*dim_mode+1:k*dim_mode ) = energy%value(j,k)%order_op(1)%Op_loc(TB_pos_start:TB_pos_end,TB_pos_start:TB_pos_end) ! energy%value(j,k) is the "k" Hamiltonian (non-zero) in the "j" shell
                enddo
            enddo

        end subroutine rewrite_H_k


#ifdef CPP_LAPACK
        subroutine diagonalise_H_k(kvector_pos, dim_mode, sense, eigval)
            implicit none
            external :: ZHEEV
            integer, intent(in) :: kvector_pos, dim_mode
            real(kind=8), intent(in) :: sense
            real(kind=8), intent(inout) :: eigval(:)

            ! Internal variable
            integer :: i

            integer :: N, INFO, N_val
            complex(kind=8), allocatable :: WORK(:), H_complex(:,:)
            real(kind=8), allocatable :: RWORK(:)

            N = size(all_E_k, 1)

            allocate( H_complex(N,N), RWORK(max(1,3*N-2)), WORK(2*N))
            ! Before diagonalising the Hamiltonian, we first have to Fourier transform it
            H_complex=get_FFT(all_E_k, all_positions, kvector_pos, dim_mode, sense)
            RWORK=0.0d0
            WORK=0.0d0
            ! diagonalising the Hamiltonian
!  Purpose
!  =======
!
!  ZHEEV computes all eigenvalues and, optionally, eigenvectors of a
!  complex Hermitian matrix A.
!
!  Arguments
!  =========
!
!  JOBZ    (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;
!          = 'V':  Compute eigenvalues and eigenvectors.
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) COMPLEX*16 array, dimension (LDA, N)
!          On entry, the Hermitian matrix A.  If UPLO = 'U', the
!          leading N-by-N upper triangular part of A contains the
!          upper triangular part of the matrix A.  If UPLO = 'L',
!          the leading N-by-N lower triangular part of A contains
!          the lower triangular part of the matrix A.
!          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!          orthonormal eigenvectors of the matrix A.
!          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!          or the upper triangle (if UPLO='U') of A, including the
!          diagonal, is destroyed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  W       (output) DOUBLE PRECISION array, dimension (N)
!          If INFO = 0, the eigenvalues in ascending order.
!
!  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The length of the array WORK.  LWORK >= max(1,2*N-1).
!          For optimal efficiency, LWORK >= (NB+1)*N,
!          where NB is the blocksize for ZHETRD returned by ILAENV.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  RWORK   (workspace) DOUBLE PRECISION array, dimension (max(1, 3*N-2))
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the algorithm failed to converge; i
!                off-diagonal elements of an intermediate tridiagonal
!                form did not converge to zero.
!
!  =====================================================================
            call ZHEEV( 'N', 'U', N, H_complex, N, eigval, WORK, 2*N, RWORK, INFO )

        end subroutine diagonalise_H_k
#endif

#ifdef CPP_INTERNAL
        subroutine diagonalise_H_k(kvector_pos, dim_mode, sense, eigval)
            use m_eigen_val_vec
            use m_invert
            implicit none
            integer :: kvector_pos, dim_mode
            real(kind=8) :: sense
            complex(kind=8), intent(out) :: eigval(:)

            ! Internal variable
            integer :: i,N

            complex(kind=8), allocatable ::  H_complex(:,:), U(:,:)

            N = size(all_E_k, 1)
            allocate(H_complex(N,N),U(N,N))
            U=0.0d0

            ! Before diagonalising the Hamiltonian, we first have to Fourier transform it
            H_complex=get_FFT(all_E_k, all_positions, kvector_pos, dim_mode, sense)
            
            call Jacobi(1.0d-7,N,H_complex,N,eigval,U,N,0)

        end subroutine diagonalise_H_k
#endif


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
            complex(kind=8), intent(in) :: E(:)
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

!        Function computing the Fermi energy of the system
        subroutine compute_Fermi_level(eps_nk, N_electrons, fermi_level, kt)
            use m_sort
            implicit none
            real(kind=8), intent(out) :: fermi_level
            real(kind=8), intent(in) :: kt
            real(kind=8), intent(in) :: eps_nk(:)
            real(kind=8), intent(in) :: N_electrons

            ! Internal variable
            integer :: i,size_ham,j
            integer, allocatable :: indices(:)
            real(kind=8) :: precision_Ef, tmp_sum

            size_ham=size(eps_nk)
            allocate( indices(size_ham) )
            call sort(size_ham, eps_nk, indices, 1.0d-5)

            write(6,'(/a,2x,2(f12.8,2x)/)') 'higher and lower eigenval',eps_nk(1),eps_nk(size_ham)

            tmp_sum = 0.0d0
            i=0
            do while ( i .le. size_ham )
              i=i+1
              fermi_level = eps_nk(i)
              tmp_sum = 0.0d0
              do j=1,i
                tmp_sum = tmp_sum + fermi_distrib(fermi_level, eps_nk(j), kt)
              enddo
              if (real(N_electrons) .le. tmp_sum) exit
            enddo
            fermi_level = eps_nk(i)

            write(6,'(/a,2x,f12.6,2x,a/)') 'fermi_level = ',  fermi_level, '[eV]'


        end subroutine compute_Fermi_level






        ! Function implementing the FD distribution for a given energy
        real(kind=8) function fermi_distrib(E_F, energy, kt)
            implicit none
            real(kind=8) :: E_F
            real(kind=8) :: kt
            real(kind=8) :: energy

            fermi_distrib=0.0d0
            if (kt.lt.1.0d-8) then
                if (energy.le.E_F) fermi_distrib=1.0d0
            else
                fermi_distrib = 1.0d0/( 1.0d0 + exp((E_F-energy)/kt ) )
            endif
        end function fermi_distrib

end module m_energy_k
