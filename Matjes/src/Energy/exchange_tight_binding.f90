module m_exchange_TB
    use m_symmetry_operators
    use m_Hamiltonian_variables, only : coeff_ham_inter_spec

    type(coeff_ham_inter_spec), target, public, protected :: exc_ham_TB

    private
    public :: get_exc_ham_TB

    contains
        subroutine get_exc_ham_TB(fname)
            use m_io_utils
            use m_io_files_utils
            use m_convert
            implicit none
            character(len=*), intent(in) :: fname
            ! internal
            integer :: io_param,neighbor_hoping
            integer :: i
            real(kind=8), allocatable :: t_local(:)

            ! Multiplicative coefficient
            exc_ham_TB%c_ham=1.0d0
            ! Name of the model
            exc_ham_TB%name='Tight-binding-Exchange'
            exc_ham_TB%N_shell=-1
            ! Order of the Hamiltonian (tensor of rank "order")
            exc_ham_TB%order=2

            io_param=open_file_read(fname)

            ! Go in the input file, search for a variable called
            ! c_tij and put its value in the variable exc_ham_TB%c_ham
            call get_parameter(io_param,fname,'c_tij',exc_ham_TB%c_ham)

            ! Count the number of hopping parameters
            call get_parameter(io_param,fname,'N_hoping',exc_ham_TB%N_shell)
            if (exc_ham_TB%N_shell.eq.-1) then
                neighbor_hoping=count_variables(io_param,'t_',fname)
            else
                neighbor_hoping=exc_ham_TB%N_shell
            endif


            ! If there are hopping parameters, read the values
            ! and put them in the allocated t_local array
            if (neighbor_hoping.ne.0) then
                allocate(t_local(neighbor_hoping))
                t_local=0.0d0

                call get_coeff(io_param,fname,'t_',t_local)
                neighbor_hoping=number_nonzero_coeff(t_local,'TB')
            endif

            if (neighbor_hoping.ne.0) exc_ham_TB%i_exist=.true.

            ! Allocate the muber of shell in exc_ham_TB
            allocate(exc_ham_TB%ham(neighbor_hoping))
            do i=1,neighbor_hoping
                allocate(exc_ham_TB%ham(i)%H(2,2))
                exc_ham_TB%ham(i)%H=0.0d0
            enddo

            ! Put the hopping parameters out of the diagonal
            if (neighbor_hoping.ne.0) then
                do i=1,neighbor_hoping
                    exc_ham_TB%ham(i)%H(1,2)=exc_ham_TB%c_ham*t_local(i)
                    exc_ham_TB%ham(i)%H(2,1)=exc_ham_TB%c_ham*t_local(i)
                enddo
            endif
        end subroutine get_exc_ham_TB
end module m_exchange_TB
