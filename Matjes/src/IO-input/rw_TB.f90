module m_rw_TB
    use m_io_utils
    use m_io_files_utils

    type parameters_TB
        real(kind=8), allocatable :: hopping(:,:)
        real(kind=8), allocatable :: onsite(:,:)
        logical :: is_magnetic=.false.
        integer :: nb_shell=-1
        integer :: nb_orbitals=-1
        integer :: nb_spin=1
        real(kind=8) :: Jsd=0.0d0

        real(kind=8) :: N_electrons  !total number of electrons 
        real(kind=8) :: kt     !smearing of fermi-dirac distribution for fermi-energy
    end type

    type(parameters_TB), public, protected :: TB_params

    private
    public :: rw_TB, check_activate_TB, get_nb_orbitals 

    contains
        subroutine rw_TB(fname)
            implicit none
            character(len=*), intent(in) :: fname
            
            ! Internal variables
            integer :: io_input
            
            ! Read in the input file the number of t_up_1, t_down_1,
            ! t_up_2, t_down_2, t_up_3, t_down_3, etc.
            ! Store that amount in nb_t_tot_spin
            ! ===> nb_t_tot_spin/2 is the number of shells
            integer :: nb_t_up,nb_t_down,nb_t,nb_t_tot_spin

            ! Read in the input file the number of mu_up_1, mu_down_1,
            ! mu_up_2, mu_down_2, t_up_3, t_down_3, etc.
            ! Store that amount in nb_mu_tot_spin
            ! ===> nb_mu_tot_spin/2 is the number of orbitals (2 electrons per orbital)
            integer :: nb_mu_down,nb_mu_up,nb_mu,nb_mu_tot_spin
            real(kind=8) :: Jsd
            logical :: activate_mag_TB

            io_input=open_file_read(fname)


            TB_params%N_electrons=1.0d0 !might want to change default
            call get_parameter(io_input, 'input', 'N_electrons', TB_params%N_electrons)
            TB_params%kt=1.0d-3 !might want to change default
            call get_parameter(io_input, 'input', 'fermi_kt',TB_params%kt)

            ! Check if magnetism is activated
            activate_mag_TB=.false.
            call get_parameter(io_input,fname,'activate_mag_TB',activate_mag_TB)

            ! Count the occurrences of different parameters
            nb_t_up=0
            nb_t_down=0
            nb_t_tot_spin=0
            nb_t=0
            if (activate_mag_TB) then
                nb_t_up=count_variables(io_input,'t_up_',fname)
                nb_t_down=count_variables(io_input,'t_down_',fname)
                nb_t_tot_spin=nb_t_up+nb_t_down
            else
                nb_t=count_variables(io_input,'t_',fname)
            endif

            nb_mu_up=0
            nb_mu_down=0
            nb_mu_tot_spin=0
            nb_mu=0
            if (activate_mag_TB) then
                nb_mu_up=count_variables(io_input,'mu_up_',fname)
                nb_mu_down=count_variables(io_input,'mu_down_',fname)
                nb_mu_tot_spin=nb_mu_up+nb_mu_down
            else
                nb_mu=count_variables(io_input,'mu_',fname)
            endif

            ! If magnetism IS activated but count(t_up) .neq. count(t_down),
            ! is is wrong ===> stop the execution
            if (activate_mag_TB .and. (nb_t_up.ne.nb_t_down)) return
            
            ! If magnetism IS activated but there are no mu_up and mu_down,
            ! is is wrong ===> stop the execution
            if (activate_mag_TB .and. (nb_mu_up.ne.nb_mu_down)) return
            
            ! If magnetism IS NOT activated but there are some t_up and/or t_down,
            ! is is wrong ===> stop the execution
            if ((.not.activate_mag_TB) .and. (nb_t_tot_spin.ne.0)) return
            
            ! If magnetism IS NOT activated but there are some mu_up and/or mu_down,
            ! is is wrong ===> stop the execution
            if ((.not.activate_mag_TB) .and. (nb_mu_tot_spin.ne.0)) return

            if (activate_mag_TB) then
                allocate(TB_params%hopping(2,nb_t_tot_spin/2))
                allocate(TB_params%onsite(2,nb_mu_tot_spin/2))
                TB_params%is_magnetic = .true.
                TB_params%nb_shell = nb_t_tot_spin/2
                TB_params%nb_orbitals = nb_mu_tot_spin/2
                TB_params%nb_spin = 2
            else
                allocate(TB_params%hopping(1,nb_t))
                allocate(TB_params%onsite(1,nb_mu))
                TB_params%is_magnetic = .false.
                TB_params%nb_shell = nb_t
                TB_params%nb_orbitals = nb_mu
                TB_params%nb_spin = 1
            endif
            TB_params%hopping=0.0d0
            TB_params%onsite=0.0d0


            ! Storage:
            !    1st line ===> 'up' values
            !    2nd line ===> 'down' values
            if (nb_t_tot_spin.ne.0) then
                call get_coeff(io_input,fname,'t_up_',TB_params%hopping(1,:))
                call get_coeff(io_input,fname,'t_down_',TB_params%hopping(2,:))
            else
                call get_coeff(io_input,fname,'t_',TB_params%hopping(1,:))
            endif
            
            if (nb_mu_tot_spin.ne.0) then
                call get_coeff(io_input,fname,'mu_up_',TB_params%onsite(1,:))
                call get_coeff(io_input,fname,'mu_down_',TB_params%onsite(2,:))
            else
                call get_coeff(io_input,fname,'mu_',TB_params%onsite(1,:))
            endif

            call close_file(fname,io_input)
        end subroutine
        
        ! Function checking if TB is active or not
        logical function check_activate_TB()
            implicit none
            check_activate_TB=.false.
            if ((TB_params%nb_shell.gt.0) .or. (TB_params%nb_orbitals.gt.0)) check_activate_TB=.true.
        end function
        
        integer function get_nb_orbitals()
            implicit none
            get_nb_orbitals=0
            if (TB_params%is_magnetic) then
                get_nb_orbitals=2*TB_params%nb_orbitals
            else
                get_nb_orbitals=TB_params%nb_orbitals
            endif
        end function
end module  m_rw_TB
