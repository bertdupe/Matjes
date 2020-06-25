module m_energy_k
    use m_energy_commons, only : energy
    use m_fftw
    implicit none

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
        subroutine compute_Fermi_level(eps_nk, N_electrons, E_f, kt)
            !subroutine which calculates the Fermi energy
            !first guess is zero temperature, difference between lowest unoccupied and higher occupied energy
            !then iterate to achieve the aimed occupation using the fermi-dirac distribution
            use m_sort
            implicit none
            real(kind=8), intent(out) :: E_f
            real(kind=8), intent(in) :: kt          !k_b*T for fermi dirac distribution
            real(kind=8), intent(in) :: eps_nk(:,:) !assume [Ns,Nk] with Nk=number of k-points, Ns=number of states
            real(kind=8), intent(in) :: N_electrons !number of states to be occupied at each k

            ! Internal variable
            integer     ::  Nk,Ns
            integer     :: i
            integer, allocatable :: indices(:)
            real(kind=8),allocatable :: tmp_E(:) !sorted energy array

            !fermi-dirac temporary variables
            !!variables determine considered energy range
            real(8),parameter     ::  cutoff_kt=10.0d0 !how many instances of kt are included in the considered energy range
                                      !arbitraily choose 10 times kt, which should be large enough n_f(-10)=0.99995
            real(8)               ::  min_E,max_E !minimal/maximal energy for calculating the fermi energy
            integer               ::  min_ind,max_ind
            !!variables within occupation convergence loop
            integer     ::  weight_ind  
            real(8)     ::  weight_aim !fermi dirac occupation aim ( in numbers of occupied states in tmp_E(min_ind:max_ind)
            real(8)     ::  dE
            !!parameters for convergence loop control
            integer,parameter     ::  fd_loop_max=50 !maximal number of fermi-dirac iterations
            real(8),parameter     ::  occ_cutoff=1.0d-4 !accurancy of the fermi dirac occupation  
            !components for first derivative:
            real(8),parameter     ::  c_deriv(7)=[1.0d0,-9.0d0,45.0d0,0.0d0,-45.0d0,9.0d0,-1.0d0]/60.0d0 
            !temporary iteration parameters for obtaining the fermi energy from the fermi-dirac distribution
            real(8)    :: E_F_ext(2),weight_ext(2)
            real(8)    :: E_F_tmp,weight_tmp

            write(*,'(A)') "Start compute fermi level"
            Ns=size(eps_nk(:,1))
            Nk=size(eps_nk(1,:))

            !get sorted eigenvalues in tmp_E
            allocate(tmp_E(Nk*Ns))
            tmp_E=reshape(eps_nk,[Nk*Ns])
            allocate(indices(Nk*Ns))
            call sort(Ns*Nk, tmp_E, indices, 1.0d-5)

            !trivial guess for fermi energy
            if(Nk*N_electrons+1 > size(tmp_E)) STOP 'Too many electrons to calculate fermi energy, reduce N_electrons?'
            E_f=(tmp_E(Nk*N_electrons)+tmp_E(Nk*N_electrons+1))*0.5d0

            !Use Fermi-dirac distribution
            !!prepare considered energy range
            !!will assume all states below min_ind are fully occupied and above max_ind are not occupied
            min_E=E_F-cutoff_kt*kt
            max_E=E_F+cutoff_kt*kt
            min_ind=1
            do i=Nk*N_electrons+1,1,-1
                if(tmp_E(i) <= min_E)then !could be done faster with an temporary array
                    min_ind=i
                    exit
                endif
            enddo
            max_ind=Nk*Ns
            do i=Nk*N_electrons,Nk*Ns
                if(tmp_E(i) >= max_E)then
                    max_ind=i
                    exit
                endif
            enddo
            weight_aim=real(Nk*N_electrons-min_ind+1,kind=8)

            !Get initial fermi energy guesses and weights 
            weight_ind=floor(weight_aim)
            if (weight_ind-3<1 .or. weight_ind+3> Nk*Ns)then
                write(*,'(A)') "Warning, very few energy values above or below the aimed occupied state"
                write(*,'(A)') "Fermi-Dirac implementation not sensible, will use initial guess Fermi energy"
                write(6,'(/a,2x,f12.6,2x,a/)') 'fermi_level = ',  E_f, '[eV]'
                return
            endif
            dE=max(dot_product(c_deriv,tmp_E(weight_ind-3:weight_ind+3)),1.0e-2)
            E_F_ext=[E_F-dE,E_F+dE]
            weight_ext(1)=fermi_weight_sum(E_F_ext(1),tmp_E(min_ind:max_ind),kt)-weight_aim
            weight_ext(2)=fermi_weight_sum(E_F_ext(2),tmp_E(min_ind:max_ind),kt)-weight_aim
            weight_ext=weight_ext

            !make sure the E_F_ext correspond to weights below and above the aim
            do while (weight_ext(1)>0.0d0 .or. weight_ext(2)<0.0d0)
                if(weight_ext(1)>0.0d0)then
                    weight_tmp=weight_ext(1)
                    E_F_tmp=E_F_ext(1)
                    E_F_ext(1)=E_F_ext(1)-(E_F_ext(2)-E_F_ext(1))*(weight_ext(1)+2.0d0)
                    weight_ext(1)=fermi_weight_sum(E_F_ext(1),tmp_E(min_ind:max_ind),kt)-weight_aim
                    weight_ext(2)=weight_tmp
                    E_F_ext(2)=E_F_tmp
                endif
                if(weight_ext(2)<0.0d0)then
                    weight_tmp=weight_ext(2)
                    E_F_tmp=E_F_ext(2)
                    E_F_ext(2)=E_F_ext(2)+(E_F_ext(2)-E_F_ext(1))*(-weight_ext(2)+2.0d0)
                    weight_ext(2)=fermi_weight_sum(E_F_ext(2),tmp_E(min_ind:max_ind),kt)-weight_aim
                    weight_ext(1)=weight_tmp
                    E_F_ext(1)=E_F_tmp
                endif
            end do

            !sucessively divide E_F and update E_F_ext and the weights until weight threshold is reached
            i=1
            weight_tmp=occ_cutoff+1.0d0
            do while (abs(weight_tmp)>occ_cutoff)
                E_F_tmp=sum(E_F_ext)*0.5d0
                weight_tmp=fermi_weight_sum(E_F_tmp,tmp_E(min_ind:max_ind),kt)-weight_aim
                if(weight_tmp<0.0d0)then
                    weight_ext(1)=weight_tmp
                    E_F_ext(1)=E_F_tmp
                else
                    weight_ext(2)=weight_tmp
                    E_F_ext(2)=E_F_tmp
                endif
                i=i+1
                if(i>fd_loop_max)then
                    write(*,'(A)') "warning, fermi dirac occupation has not reached the wanted cutoff"
                    write(*,'(A,E16.8)') "weight-difference",weight_tmp
                    write(*,'(A)') "Possibly continuing with wrong Fermi energy"
                    exit
                endif
            end do 
            E_f=E_F_tmp
            write(6,'(/a,2x,f12.6,2x,a/)') 'fermi_level = ',  E_f, '[eV]'

        end subroutine compute_Fermi_level


        function fermi_weight_sum(E_F,energy,kt)result(weight)
            real(8)     :: weight
            real(8)     :: E_F
            real(8)     :: kt
            real(8)     :: energy(:)

            integer     ::  i
            weight=0.0d0
            do i=1,size(energy)
                weight=weight+fermi_distrib(E_F, energy(i), kt)
            end do
        end function



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
                fermi_distrib = 1.0d0/( 1.0d0 + exp((energy-E_F)/kt ) )
            endif
        end function fermi_distrib

end module m_energy_k
