module m_bandstructure

    use m_basic_types, only : vec_point
    use m_energy_commons, only : energy
    use m_fftw

    ! all vectors on one line (must be updated at each line)
    complex(kind=16), allocatable, dimension(:) :: all_vectors ! will contain all the order parameters in one line
    real(kind=8), allocatable :: all_positions(:,:) ! array containing the coordinates of all r-r'

    ! all E matrix on one line (must be done once since the table is always done in the same order)
    ! The ordering is as follows: line "i" corresponds to site "i". Each column of a given line
    ! corresponds to the different Hamiltonians: column "i" is the onsite, column "i+1" is the
    ! first shell Hamiltonian, column "i+2" is the 2nd shell Hamiltonian, etc.
    complex(kind=16), allocatable, dimension(:,:) :: all_E


    interface print_band_struct
      module procedure print_band_struct_1e,print_band_struct_Ne
    end interface print_band_struct
    private
    public :: set_E_bandstructure,calculate_dispersion,print_band_struct

    contains
        subroutine set_E_bandstructure(dim_mode,pos)
            implicit none
            integer, intent(in) :: dim_mode
            real(kind=8), intent(in) :: pos(:,:)

            ! Internal variables
            integer :: N,i,i_vois

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
                        all_vectors((k-1)*dim_mode+1:k*dim_mode) = cmplx(all_mode(i_vois)%w, kind=16) !all_mode is the vector of ordre parameter for each site. We put everything in a big vector
                    enddo
                    dispersion(i) = dispersion(i) + get_FFT(all_E, all_positions, all_vectors, cmplx(all_mode(j)%w, kind=16), i, dim_mode*N, dim_mode, -1.0d0) !dispersion(i)=dispersion(i)+... accounts for the sum over r
                enddo
            enddo
            dispersion = dispersion/real(N_cell)
        end subroutine calculate_dispersion

        ! print band structure
        subroutine print_band_struct_1e(fname,dispersion)
        use m_io_utils
        use m_io_files_utils
        implicit none
        complex(kind=16), intent(inout) :: dispersion(:)
        character(len=*), intent(in) :: fname
        ! internal
        integer :: io_band,i,N_k

        N_k=size(dispersion)

        io_band=open_file_write(fname)

        do i=1,N_k
           write(io_band,'(2(E20.12E3,3x))') real(dispersion(i)),aimag(dispersion(i))
        enddo
        call close_file(fname,io_band)

        end subroutine


        ! print band structure
        subroutine print_band_struct_Ne(fname,dispersion)
        use m_io_utils
        use m_io_files_utils
        use m_convert
        implicit none
        real(kind=8), intent(inout) :: dispersion(:,:)
        character(len=*), intent(in) :: fname
        ! internal
        integer :: io_band,i,N_k,j,N_e
        character(len=30) :: form

        N_k=size(dispersion,2)
        N_e=size(dispersion,1)
        form=convert('(',N_e,'(E20.12E3,3x))')

        io_band=open_file_write(fname)

        do i=1,N_k
           write(io_band,form) (dispersion(j,i),j=1,N_e)
        enddo
        call close_file(fname,io_band)

        end subroutine

end module m_bandstructure
