module m_chem_pot_TB
    use m_symmetry_operators
    use m_Hamiltonian_variables, only : coeff_ham_inter_spec
    use m_rw_TB
    use m_lattice, only : my_order_parameters

    type(coeff_ham_inter_spec), target, public, protected :: onsite_ham_TB

    private
    public :: get_onsite_ham_TB

    contains
        subroutine get_onsite_ham_TB(fname,dim_ham)
            use m_io_utils
            use m_io_files_utils
            use m_convert
            implicit none
            character(len=*), intent(in) :: fname
            integer, intent(in) :: dim_ham
            
            ! Internal variables
            integer :: io_param,chem_pot_count
            character(len=50) :: form
            integer :: i, j, x_start, x_end
            real(kind=8), allocatable :: chem_pot_local(:)
            
            do i=1, size(my_order_parameters)
                if (my_order_parameters(i)%name.eq.'Tight-binding') then
                    x_start=my_order_parameters(i)%start
                    x_end=my_order_parameters(i)%end
                endif
            enddo
            
            chem_pot_count=TB_params%nb_orbitals

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

            if (chem_pot_count.ne.0) then
              onsite_ham_TB%i_exist=.true.
            else
              return
            endif

            ! Allocate the different blocs in the total Hamiltonian
            allocate(onsite_ham_TB%ham(1))
            allocate(onsite_ham_TB%ham(1)%H(dim_ham,dim_ham))
            onsite_ham_TB%ham(1)%H=0.0d0

            allocate(chem_pot_local(x_end-x_start+1))
            chem_pot_local=reshape( TB_params%onsite, (/x_end-x_start+1/) )
            j=0           
            do i=x_start,x_end
                j=j+1
                onsite_ham_TB%ham(1)%H(i,i)=onsite_ham_TB%c_ham*chem_pot_local(j)
            enddo
            
            form=convert('(',dim_ham,'(f12.8,2x))')
            write(6,'(a)') ''
            write(6,'(a)') 'Onsite tight-binding is OK'
            do i=1,dim_ham
                write(6,form) onsite_ham_TB%ham(1)%H(:,i)
            enddo
            write(6,'(a)') ''
        end subroutine get_onsite_ham_TB
end module m_chem_pot_TB
