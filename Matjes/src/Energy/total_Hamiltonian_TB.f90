module m_total_Hamiltonian_TB
    use m_exchange_TB
    use m_chem_pot_TB
    use m_Hamiltonian_variables, only: coeff_ham_inter_spec
    use m_convert
    use m_J_sd_exchange

    type(coeff_ham_inter_spec), target, public, protected :: ham_tot_TB

    private
    public :: get_total_Hamiltonian_TB

    contains
        subroutine get_total_Hamiltonian_TB(fname,dim_ham)
            implicit none
            character(len=*), intent(in) :: fname
            integer, intent(in) :: dim_ham

            ! Internal variables
            integer :: size_ham
            integer :: i,j
            character(len=50) :: form


            call get_onsite_ham_TB(fname,dim_ham)
            call get_exc_ham_TB(fname,dim_ham)
            call get_Jsd_ham_TB(fname,dim_ham)

            if (.not. (onsite_ham_TB%i_exist .or. exc_ham_TB%i_exist) ) return

            ham_tot_TB%c_ham=1.0d0
            ham_tot_TB%name='Total-Tight-Binding'
            ham_tot_TB%i_exist=.true.
            ham_tot_TB%order=2

            ! size_ham is the size of the total Hamiltonian. We need to
            ! now the number of shells that are considered, which is realised
            ! thanks to size(exc_ham_TB%ham).
            ! We add 1 because the onsite Hamiltonian has only 1 shell
            size_ham=size(exc_ham_TB%ham)+1
            ham_tot_TB%N_shell=size_ham
            allocate(ham_tot_TB%ham(size_ham))

            do i=1,size_ham
                allocate(ham_tot_TB%ham(i)%H(dim_ham,dim_ham))
                ham_tot_TB%ham(i)%H=0.0d0
            enddo

            ! Filling of total Hamiltonian matrix
            if (onsite_ham_TB%i_exist) ham_tot_TB%ham(1)%H = ham_tot_TB%ham(1)%H + onsite_ham_TB%ham(1)%H
            if (exc_ham_TB%i_exist) then
              do i=2,size_ham
                 ham_tot_TB%ham(i)%H = ham_tot_TB%ham(i)%H + exc_ham_TB%ham(i-1)%H
              enddo
            endif
            !if (Jsd_ham_TB%i_exist) ham_tot_TB%ham(1)%H = ham_tot_TB%ham(1)%H + Jsd_ham_TB%ham(1)%H
            
            form=convert('(',dim_ham,'(f12.8,2x))')
            write(6,'(a)') ''
            write(6,'(a)') 'Total Hamiltonian tight-binding is OK'
            do j=1,size(ham_tot_TB%ham)
                write(6,'(a,I8)') 'Now in shell ', j-1
                do i=1,dim_ham
                    write(6,form) ham_tot_TB%ham(j)%H(:,i)
                enddo
            enddo
            write(6,'(a)') ''
        end subroutine get_total_Hamiltonian_TB
end module m_total_Hamiltonian_TB
