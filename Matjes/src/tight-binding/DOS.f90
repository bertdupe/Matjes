module m_DOS
    ! Used to read an input file
    use m_io_utils
    use m_io_files_utils

    real(kind=8) :: E_F = 0.0d0
    real(kind=8) :: from = -1.0d0
    real(kind=8) :: to = 1.0d0
    real(kind=8) :: smearing = 1.0d-2
    integer :: n_pt = 1000
    character(len=20) :: smearing_type = 'gaussian'

    complex(kind=16), allocatable :: DOS(:),E_DOS(:)

    private
    public :: read_params_DOS, init_Evector_DOS, compute_DOS, print_DOS
    contains
        subroutine read_params_DOS(fname)
            implicit none
            character(len=*), intent(in) :: fname

            ! Internal variables
            integer :: io_input

            io_input=open_file_read(fname)

            call get_parameter(io_input, fname, 'E_F', E_F)
            call get_parameter(io_input, fname, 'n_pt', n_pt)
            call get_parameter(io_input, fname, 'from', from)
            call get_parameter(io_input, fname, 'to', to)
            call get_parameter(io_input, fname, 'smearing_type', smearing_type)
            call get_parameter(io_input, fname, 'smearing', smearing)

            call close_file(fname, io_input)
        end subroutine read_params_DOS



        ! Subroutine initialising the energy vector required by the DOS
        subroutine init_Evector_DOS()
            implicit none

            ! Internal variable
            integer :: i
            real(kind=8) :: range_values

            allocate( E_DOS(n_pt) )
            E_DOS = 0.0d0

            range_values = to - from
            do i=1, n_pt
                E_DOS(i) = cmplx(from, kind=16) + cmplx( (i-1)*range_values/(n_pt - 1), kind=16)
            enddo
        end subroutine init_Evector_DOS



        ! Function computing the DOS
        subroutine compute_DOS(disp_en, N_cell)
            implicit none
            integer, intent(in) :: N_cell
            complex(kind=16), intent(in) :: disp_en(:)

            ! Internal variable
            integer :: i, j
            real(kind=8) :: lower_bound,upper_bound

            allocate( DOS(n_pt) )
            DOS = 0.0d0

            do i = 1, size(disp_en)
                do j = 1, size( E_DOS )
                    lower_bound=0.7*real(disp_en(i))*exp( -real(E_DOS(j)-disp_en(i))**2/(2*smearing**2))
                    upper_bound=1.3*real(disp_en(i))*exp( -real(E_DOS(j)-disp_en(i))**2/(2*smearing**2))
                    if( ((real(E_DOS(j)) .le. upper_bound)).and.((real(E_DOS(j)) .ge. lower_bound)) )  DOS(j)=DOS(j)+1.0d0
                enddo
            enddo
            DOS = DOS/real(N_cell)
        end subroutine compute_DOS

        ! Function printing the DOS
        subroutine print_DOS(fname)
        use m_io_utils
        use m_io_files_utils
        implicit none
        character(len=*), intent(in) :: fname
        ! internal
        integer :: io_DOS,i,N_k

        io_DOS=open_file_write(fname)
        do i=1,n_pt
           write(io_DOS,'(3(E20.12E3,3x))') real(E_DOS(i)),real(DOS(i)),aimag(DOS(i))
        enddo
        call close_file(fname,io_DOS)

        end subroutine

end module m_DOS
