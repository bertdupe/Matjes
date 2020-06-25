module m_wavefunction
    use m_basic_types, only : vec_point

    private
    public :: check_norm_wavefct
    contains
        ! Function computing the normalisation of the wavefunction. It performs
        !   sum_i c_i^* X c_i
        ! and should return the number of particles.
        ! Input:
        !   _ all_mode: vector containing the ordre parameter for each site
        !   _ dim_mode: the length of the ordre parameter for each site
        ! Output:
        !   _ N_electrons: total number of electrons in the system
        subroutine check_norm_wavefct(all_mode, N_electrons)
            implicit none
            real(kind=8), intent(out) :: N_electrons
            type(vec_point), intent(in) :: all_mode(:)

            ! Internal variables
            integer :: i

            N_electrons = 0.0d0
            do i = 1, size(all_mode)
                N_electrons = N_electrons + sum(all_mode(i)%w**2)
                write(*,*) i,N_electrons
            enddo
        end subroutine check_norm_wavefct

end module m_wavefunction
