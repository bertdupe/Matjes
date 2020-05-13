module m_chem_pot_TB
    use m_symmetry_operators
    use m_Hamiltonian_variables, only : coeff_ham_inter_spec

    type(coeff_ham_inter_spec), target, public, protected :: onsite_ham_TB

    private
    public :: get_onsite_ham_TB

    contains
        subroutine get_onsite_ham_TB(fname)
            use m_io_utils
            use m_io_files_utils
            use m_convert
            implicit none
            character(len=*), intent(in) :: fname
            integer :: io_param,chem_pot_count
            real(kind=8), allocatable :: chemPot_local(:)
            integer :: i

            ! Multiplicative coefficient
            onsite_ham_TB%c_ham=1.0d0
            ! Name of the model
            onsite_ham_TB%name='Tight-binding-Onsite'
            onsite_ham_TB%N_shell=-1
            ! Order of the Hamiltonian (tensor of rank "order")
            onsite_ham_TB%order=2

            io_param=open_file_read(fname)

            ! Go in the input file, search for a variable called
            ! c_Ei and put its value in the variable onsite_ham_TB%c_ham
            call get_parameter(io_param,fname,'c_Ei',onsite_ham_TB%c_ham)

            ! Count the number of different chemical potentials
            chem_pot_count=count_variables(io_param,'mu_',fname)

            ! If there are muliple values for the chemical potential,
            ! put each of these values in a vector
            if (chem_pot_count.ne.0) then
                allocate(chemPot_local(chem_pot_count))
                chemPot_local=0.0d0

                call get_coeff(io_param,fname,'mu_',chemPot_local)
                chem_pot_count=number_nonzero_coeff(chemPot_local,'onsite_ham_TB')
            endif

            if (chem_pot_count.ne.0) onsite_ham_TB%i_exist=.true.

            ! Allocate the different blocs in the total Hamiltonian
            allocate(onsite_ham_TB%ham(chem_pot_count))
            do i=1,chem_pot_count
                allocate(onsite_ham_TB%ham(i)%H(2,2))
                onsite_ham_TB%ham(i)%H=0.0d0
            enddo

            if (chem_pot_count.ne.0) then
                do i=1,chem_pot_count
                    onsite_ham_TB%ham(i)%H(1,1)=onsite_ham_TB%c_ham*chemPot_local(i)
                    onsite_ham_TB%ham(i)%H(2,2)=onsite_ham_TB%c_ham*chemPot_local(i)
                    !write(*,*) TB%ham(i)%H
                enddo
            endif
        end subroutine get_onsite_ham_TB
end module m_chem_pot_TB
