module m_wavefunction
    use m_basic_types, only : vec_point

    private
    public :: check_norm_wavefct, compute_Fermi_energy
    contains
        ! Function computing the normalisation of the wavefunction. It performs
        !   sum_i c_i^* X c_i
        ! and should return the number of particles.
        ! Input:
        !   _ all_mode: vector containing the ordre parameter for each site
        !   _ dim_mode: the length of the ordre parameter for each site
        subroutine check_norm_wavefct(all_mode, dim_mode)
            implicit none
            integer, intent(in) :: dim_mode
            type(vec_point), intent(in) :: all_mode(:)

            ! Internal variables
            integer :: i
            real(kind=8) :: nb_particles

            nb_particles = 0.0d0
            do i = 1, size(all_mode)
                nb_particles = nb_particles + sum(all_mode(i)%w(4:dim_mode))*sum(all_mode(i)%w(4:dim_mode))
            enddo
            write(*,*) 'In FILE ', __FILE__
            write(*,*) 'nb_particles = ', nb_particles
        end subroutine check_norm_wavefct



        ! Subroutine computing the Fermi energy.
        ! The Fermi energy is the highest populated energy level @ 0K.
        subroutine compute_Fermi_energy()
            implicit none
        end subroutine compute_Fermi_energy
end module m_wavefunction