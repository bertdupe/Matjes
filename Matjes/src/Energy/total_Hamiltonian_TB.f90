module m_total_Hamiltonian_TB
    use m_exchange_TB
    use m_chem_pot_TB
    use m_Hamiltonian_variables, only: coeff_ham_inter_spec

    type(coeff_ham_inter_spec), target, public, protected :: ham_tot_TB

    private
    public :: get_total_Hamiltonian

    contains
        subroutine get_total_Hamiltonian(fname)
            implicit none
            character(len=*), intent(in) :: fname

            ! Internal variables
            integer :: size_ham, shape_ham(2)
            integer :: i


            call get_onsite_ham_TB(fname)
            call get_exc_ham_TB(fname)

            if (.not. (onsite_ham_TB%i_exist .or. exc_ham_TB%i_exist) ) return

            ham_tot_TB%c_ham=1.0d0
            ham_tot_TB%name='Total-Tight-Binding'
            ham_tot_TB%order=2

            size_ham=size(exc_ham_TB%ham)+size(onsite_ham_TB%ham)
            ham_tot_TB%N_shell=size_ham
            allocate(ham_tot_TB%ham(size_ham))

            shape_ham = shape(exc_ham_TB%ham(1)%H)
            do i=1,size_ham
                allocate(ham_tot_TB%ham(i)%H(shape_ham(1), shape_ham(2)))
                ham_tot_TB%ham(i)%H=0.0d0
            enddo

            ! Filling of total Hamiltonian matrix
            do i=1,size_ham
                !if (onsite_ham_TB%i_exist .and. exc_ham_TB%i_exist) then
                if (onsite_ham_TB%i_exist .and. (i==1)) ham_tot_TB%ham(i)%H =ham_tot_TB%ham(i)%H+onsite_ham_TB%ham(i)%H
                if (exc_ham_TB%i_exist .and. (i>1)) ham_tot_TB%ham(i)%H = ham_tot_TB%ham(i)%H+exc_ham_TB%ham(i-1)%H
            enddo
        end subroutine get_total_Hamiltonian
end module m_total_Hamiltonian_TB
