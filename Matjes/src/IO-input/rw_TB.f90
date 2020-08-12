module m_rw_TB
    use m_TB_types
    use m_io_utils
    use m_io_files_utils
    implicit none

    type(parameters_TB), public,protected :: TB_params

    private
    public :: rw_TB, check_activate_TB, get_nb_orbitals ,get_parameters_io_TB

    contains

        subroutine get_parameters_io_TB(param)
            type(parameters_TB),intent(out) :: param
            param=TB_params

        end subroutine

        subroutine rw_TB(fname)
            character(len=*), intent(in) :: fname
            
            Call read_TB_H(fname,TB_params%io_H)
            Call read_TB_EF(fname,TB_params%io_ef)
            Call read_TB_dos(fname,TB_params%io_dos)
            Call read_TB_flow(fname,TB_params%flow)
            Call read_TB_highs(fname,TB_params%io_highs)

            if(TB_params%io_H%nb_shell==2) TB_params%is_mag=.True.
        end subroutine


        subroutine read_TB_highs(fname,highs)
            use m_convert
            character(len=*), intent(in)             :: fname
            type(parameters_TB_IO_highs),intent(out)     :: highs
            integer     ::  io_input

            real(8),allocatable     ::  tmp(:,:)

            io_input=open_file_read(fname)
            call get_parameter(io_input,'input','N_highsym',highs%N_highsym)
            if(highs%N_highsym < 1)then
                return
            endif
            allocate(highs%k_highs(3,highs%N_highsym),source=0.0d0)
            allocate(tmp(highs%N_highsym,3),source=0.0d0)  !necessary because the order of reading in get_parameter
            call get_parameter(io_input,'input','k_highs_pts',highs%N_highsym,3,tmp)
            highs%k_highs=transpose(tmp)
            deallocate(tmp)
            allocate(highs%k_highs_frac(highs%N_highsym),source=1)
            call get_parameter(io_input,'input','k_highs_frac',highs%N_highsym,highs%k_highs_frac)
            call get_parameter(io_input,'input','k_highs_dist',highs%aim_dist)
            call close_file('input',io_input)
        end subroutine 



        subroutine read_TB_flow(fname,flow)
            character(len=*), intent(in)             :: fname
            type(parameters_TB_IO_flow),intent(out)     :: flow
            integer     ::  io_input

            io_input=open_file_read(fname)
            call get_parameter(io_input,fname,'do_TB_r',flow%do_r)
            call get_parameter(io_input,fname,'do_dos_r',flow%dos_r)
            call get_parameter(io_input,fname,'do_occ_r',flow%occ_r)
            call get_parameter(io_input,fname,'do_spec_r',flow%spec_r)
            call get_parameter(io_input,fname,'do_fermi_r',flow%fermi_r)

            call get_parameter(io_input,fname,'do_TB_k',flow%do_k)
            call get_parameter(io_input,fname,'do_dos_k',flow%dos_k)
            call get_parameter(io_input,fname,'do_fermi_k',flow%fermi_k)
            call get_parameter(io_input,fname,'do_highs_k',flow%highs_k)
            call close_file(fname,io_input)

        end subroutine

        subroutine read_TB_dos(fname,io_dos)
            character(len=*), intent(in)             :: fname
            type(parameters_TB_IO_dos),intent(out)     :: io_dos
            integer     ::  io_input

            io_input=open_file_read(fname)
            call get_parameter(io_input,fname,'dos_sigma',io_dos%sigma)
            call get_parameter(io_input,fname,'dos_E_ext',2,io_dos%E_ext)
            call get_parameter(io_input,fname,'dos_dE',io_dos%dE)
            call close_file(fname,io_input)

        end subroutine

        subroutine read_TB_EF(fname,io_ef)
            character(len=*), intent(in)             :: fname
            type(parameters_TB_IO_EF),intent(out)     :: io_ef
            integer     ::  io_input

            io_input=open_file_read(fname)
            call get_parameter(io_input, fname, 'N_electrons', io_ef%N_electrons)
            call get_parameter(io_input, fname, 'fermi_kt', io_ef%kt)
            call get_parameter(io_input, fname, 'TB_EF', io_ef%E_F_in)
            call close_file(fname,io_input)

        end subroutine
            
        subroutine read_TB_H(fname,TB_params)
            use m_convert
            character(len=*), intent(in)             :: fname
            type(parameters_TB_IO_H),intent(out)     :: TB_params
            
            ! Internal variables
            integer :: io_input
            
            ! do loops
            integer :: i,j

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
            logical :: activate_mag_TB
            ! name of the variables
            character(len=50) :: seed_name,seed_name_mag


            io_input=open_file_read(fname)


!            TB_params%N_electrons=1.0d0 !might want to change default
!            call get_parameter(io_input, 'input', 'N_electrons', TB_params%N_electrons)
!            TB_params%kt=1.0d-3 !might want to change default
!            call get_parameter(io_input, 'input', 'fermi_kt',TB_params%kt)

            ! Check if magnetism is activated
            activate_mag_TB=.false.
            call get_parameter(io_input,fname,'activate_mag_TB',activate_mag_TB)

            ! Count the occurrences of different parameters
            nb_t_up=0
            nb_t_down=0
            nb_t_tot_spin=0
            nb_t=0
            seed_name='t_'
            if (activate_mag_TB) then
                seed_name_mag=convert(seed_name,'up_1_')
                nb_t_up=count_variables(io_input,seed_name_mag,fname)
                seed_name_mag=convert(seed_name,'down_1_')
                nb_t_down=count_variables(io_input,seed_name_mag,fname)
                nb_t_tot_spin=nb_t_up+nb_t_down
            else
                nb_t=count_variables(io_input,seed_name,fname)
            endif

            nb_mu_up=0
            nb_mu_down=0
            nb_mu_tot_spin=0
            nb_mu=0
            seed_name='mu_'
            if (activate_mag_TB) then
                seed_name_mag=convert(seed_name,'up_')
                nb_mu_up=count_variables(io_input,seed_name_mag,fname)
                seed_name_mag=convert(seed_name,'down_')
                nb_mu_down=count_variables(io_input,seed_name_mag,fname)
                nb_mu_tot_spin=nb_mu_up+nb_mu_down
            else
                nb_mu=count_variables(io_input,seed_name,fname)
            endif

            !
            ! should have at least one orbital otherwise, it stops
            !

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

            call get_parameter(io_input,fname,'nb_shell',TB_params%nb_shell)
            call get_parameter(io_input,fname,'nb_orbitals',TB_params%nb_orbitals)

            if ((TB_params%nb_shell.lt.0).or.(TB_params%nb_orbitals.lt.0)) stop 'nb_shell or nb_orbitals not read in input'

            if (activate_mag_TB) then
                allocate(TB_params%hopping(2,TB_params%nb_orbitals,TB_params%nb_shell))
                allocate(TB_params%onsite(2,TB_params%nb_orbitals))
                allocate(TB_params%Jsd(TB_params%nb_orbitals))
                TB_params%nb_spin = 2
            else
                allocate(TB_params%hopping(1,TB_params%nb_orbitals,TB_params%nb_shell))
                allocate(TB_params%onsite(1,TB_params%nb_orbitals))
                allocate(TB_params%Jsd(1))
                TB_params%nb_spin = 1
            endif
            TB_params%hopping=0.0d0
            TB_params%onsite=0.0d0
            TB_params%Jsd=0.0d0
            allocate(TB_params%delta(TB_params%nb_orbitals),source=cmplx(0.0d0,0.0d0,8))


            ! Storage:
            !    1st line ===> 'up' values
            !    2nd line ===> 'down' values
            if (nb_t_tot_spin.ne.0) then
                do j=1,2
                  if (j.eq.1) then
                    seed_name=convert('t_','up_')
                  else
                    seed_name=convert('t_','down_')
                  endif

                  do i=1,TB_params%nb_orbitals
                    seed_name_mag=convert(seed_name,i,'_')
                    call get_coeff(io_input,fname,seed_name_mag,TB_params%hopping(j,i,:))
                  enddo
                enddo
            else
                do i=1,TB_params%nb_orbitals
                  seed_name=convert('t_',i,'_')
                  call get_coeff(io_input,fname,seed_name,TB_params%hopping(1,i,:))
                enddo
            endif
            if (nb_mu_tot_spin.ne.0) then
                do j=1,2
                  if (j.eq.1) then
                    seed_name=convert('mu_','up_')
                  else
                    seed_name=convert('mu_','down_')
                  endif

                  call get_coeff(io_input,fname,seed_name,TB_params%onsite(j,:))

                enddo
            else
                call get_coeff(io_input,fname,'mu_',TB_params%onsite(1,:))
            endif

            call get_coeff(io_input,fname,'J_sd_',TB_params%Jsd)
            call get_parameter(io_input,fname,'SC_delta',TB_params%nb_orbitals,TB_params%delta)
            call get_parameter(io_input,fname,'TB_sparse',TB_params%sparse)
            call get_parameter(io_input,fname,'TB_diag',TB_params%i_diag)
            call close_file(fname,io_input)

        end subroutine
        
        ! Function checking if TB is active or not
        logical function check_activate_TB()
            implicit none
            check_activate_TB=.false.
            if ((TB_params%io_H%nb_shell.gt.0) .or. (TB_params%io_H%nb_orbitals.gt.0)) check_activate_TB=.true.
        end function
        
        integer function get_nb_orbitals()
            implicit none
            get_nb_orbitals=TB_params%io_H%nb_spin*TB_params%io_H%nb_orbitals
        end function
end module  m_rw_TB
