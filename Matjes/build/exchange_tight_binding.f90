module m_exchange_TB
    use m_symmetry_operators
    use m_Hamiltonian_variables, only : coeff_ham_inter_spec
    use m_rw_TB
    use m_lattice, only : my_order_parameters

    type(coeff_ham_inter_spec), target, public, protected :: exc_ham_TB

    private
    public :: get_exc_ham_TB

    contains
        subroutine get_exc_ham_TB(fname,dim_ham)
            use m_io_utils
            use m_io_files_utils
            use m_convert
            implicit none
            integer, intent(in) :: dim_ham
            character(len=*), intent(in) :: fname

            ! internal
            integer :: io_param,nb_shell
            character(len=50) :: form
            integer :: i, j
            integer :: x_start, x_end
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

            nb_shell=TB_params%nb_shell
            if (nb_shell.ne.0) exc_ham_TB%i_exist=.true.

            do i=1, size(my_order_parameters)
                if (my_order_parameters(i)%name.eq.'Tight-binding') then
                    x_start=my_order_parameters(i)%start
                    x_end=my_order_parameters(i)%end
                endif
            enddo

            ! Allocate the different blocs in the total Hamiltonian
            allocate(exc_ham_TB%ham(nb_shell))
            do i=1,nb_shell
                allocate(exc_ham_TB%ham(i)%H(dim_ham,dim_ham))
                exc_ham_TB%ham(i)%H=0.0d0
            enddo

            allocate(t_local(x_end-x_start+1))
            t_local=reshape( TB_params%hopping, (/x_end-x_start+1/) )

            ! Put the hopping parameters on the diagonal
            do j=1,nb_shell
                do i=x_start,x_end
                    exc_ham_TB%ham(j)%H(i,i)=exc_ham_TB%c_ham*t_local(j)
                enddo
            enddo

            form=convert('(',dim_ham,'(f12.8,2x))')
            write(6,'(a)') ''
            write(6,'(a)') 'Exchange tight-binding is OK'
            do j=1,size(exc_ham_TB%ham)
              write(6,'(a,I8)') 'Exchange tight-binding - shell',j
              do i=1,dim_ham
                write(6,form) exc_ham_TB%ham(j)%H(:,i)
              enddo
            enddo
            write(6,'(a)') ''

        end subroutine get_exc_ham_TB
end module m_exchange_TB 