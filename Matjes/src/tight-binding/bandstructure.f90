module m_bandstructure
    use m_basic_types, only : vec_point
    use m_energy_commons, only : energy
    use m_fftw

    ! all vectors on one line (must be updated at each line)
    complex(kind=16), allocatable, dimension(:) :: all_vectors
    real(kind=8), allocatable :: all_positions(:,:) ! array containing the coordinates of all r-r'

    ! all E matrix on one line (must be done once since the table is always done in the same order)
    ! The ordering is as follows: line "i" corresponds to site "i". Each column of a given line
    ! corresponds to the different Hamiltonians: column "i" is the onsite, column "i+1" is the
    ! first shell Hamiltonian, column "i+2" is the 2nd shell Hamiltonian, etc.
    complex(kind=16), allocatable, dimension(:,:) :: all_E

    public :: set_E_bandstructure,calculate_dispersion

    contains
        subroutine set_E_bandstructure(dim_mode,pos)
            implicit none
            integer, intent(in) :: dim_mode
            real(kind=8), intent(in) :: pos(:,:)

            ! Internal variables
            integer :: N,i,j,i_vois

            ! Gives the number of sites in the lattice
            N=size(energy%line(:,1))

            allocate(all_vectors(dim_mode*N),all_E(dim_mode,dim_mode*N), all_positions(3,N))
            all_vectors=cmplx(0.0d0, kind=16)
            all_E=cmplx(0.0d0, kind=16)
            all_positions=0.0d0

            do i=1,N
                i_vois=energy%line(i,1) !all neighbours of site 1
                all_E(:,(i-1)*dim_mode+1:i*dim_mode)=transpose(cmplx(energy%value(i,1)%order_op(1)%Op_loc, kind=16))
                all_positions(:,i) = pos(:,i_vois) !could have written all_positions(:,i) = pos(:,i_vois)-pos(:,1), but pos(:,1) is always 0
            enddo
        end subroutine set_E_bandstructure



        ! Function computing the totyal energy for each kpoint in the mesh
        subroutine calculate_dispersion(all_mode, dispersion, dim_mode, nb_kpoint, N_cell)
            implicit none
            type(vec_point),intent(in) :: all_mode(:)
            integer, intent(in) :: dim_mode, nb_kpoint, N_cell
            complex(kind=16), intent(inout) :: dispersion(:)

            ! Internal variable
            integer :: i, j, k, N, i_vois

            ! Gives the number of sites in the lattice
            N=size(energy%line(:,1))

            do i=1, nb_kpoint
                do j=1, N_cell !loop to get the sites
                    do k=1, N !loop to get the neighbours of site j
                        i_vois=energy%line(k,j) !contain the neighbour "k" of site "j"
                        all_vectors((k-1)*dim_mode+1:k*dim_mode)=cmplx(all_mode(i_vois)%w,kind=16) !all_mode is the vector of ordre parameter for each site. We put everything in a big vector
                    enddo
                    dispersion(i) = dispersion(i) + get_FFT(all_E, all_positions, all_vectors, cmplx(all_mode(j)%w, kind=16), i, dim_mode*N, dim_mode, -1.0d0) !dispersion(i)=dispersion(i)+... accounts for the sum over r
                enddo
            enddo
            dispersion = dispersion/real(N_cell)
        end subroutine calculate_dispersion



        subroutine initiate_input_E(fermi_energy, in_en)
            implicit none
            real(kind=8), intent(in) :: fermi_energy
            complex(kind=16), intent(inout) :: in_en(:)

            ! Internal variables
            integer :: i

            do i=1, size(in_en)
                in_en(i) = complex( (i-1)*0.75*fermi_energy/(size(in_en)-1), 0.0d0)
            enddo
        end subroutine initiate_input_E



        ! Computes the DOS according to the formula
        !   DOS(E) = (1/V) * sum_i Kronecker( E - E(k_i) )
        ! Applying strictly the definition is perhaps arsh
        !   ===> perhaps use a gaussian smearing
        ! Inputs:
        !   _ fermi_energy: maximum energy avalaible (above, DOS=0)
        !   _ disp_en: array of energy given by the computation of the dispersion relation
        !   _ in_en: array of energy values for the DOS
        ! Output:
        !   _ DOS: DOS of the system
!        subroutine compute_DOS(disp_en, in_en, DOS, N_cell)
!            implicit none
!            integer, intent(in) :: N_cell
!            complex(kind=16), intent(in) :: disp_en(:), in_en(:)
!            complex(kind=16), intent(inout) :: DOS(:)
!
!            ! Internal variables
!            integer :: i, j, nb_E_levels
!            real(kind=8) :: dE, sigma
!
!            sigma=1.0d-2
!
!            do i=1, size(disp_en)
!                do j=1, size(in_en)
!!                    if ( (real(in_en(j))-dE.le.real(disp_en(i))) .and. (real(in_en(j))+dE.ge.real(disp_en(i))) ) DOS(j)=DOS(j)+1
!                    if( ((real(in_en(j)) .le. 1.3*real(disp_en(i)*exp( -((in_en(j)-disp_en(i))**2)/(2*sigma**2)))))\
!                            .and.\
!                            ((real(in_en(j)).ge.real(0.7*disp_en(i)*exp( -((in_en(j)-disp_en(i))**2)/(2*sigma**2))))) )\
!                            DOS(j)=DOS(j)+1
!                enddo
!            enddo
!            DOS = DOS/cmplx(N_cell, kind=16)
!        end subroutine compute_DOS
end module m_bandstructure