module m_rw_H
    use m_input_H_types
    public

    contains
    subroutine rw_H(H_io)
        use m_io_files_utils, only : open_file_read,close_file
        use m_anisotropy_heisenberg, only: read_anisotropy_input
        use m_zeeman, only: read_zeeman_input
        use m_exchange_heisenberg_J, only: read_J_input
        use m_exchange_heisenberg_D, only: read_D_input
        use m_exchange_TJ, only: read_TJ_input
        use m_coupling_ME_J, only: read_ME_J_input
        use m_coupling_ME_D, only: read_ME_D_input
        use m_harmonic_phonon, only : read_F_input
        use m_ASR_phonon, only : read_ASR_Ph_input
        use m_stark, only : read_stark_input
        use m_Mag_Biq, only : read_Mag_Biq_input
        use m_4spin, only: read_sp4_input
        use m_dipolar_magnetic, only: read_dip_input
        use m_exchange_heisenberg_general, only : read_ExchG_input
        use m_spincurrent, only : read_SC_input
        use m_phonon_rank4, only : read_PH4_input
        use m_Ph_Biq, only : read_Ph_Biq_input
        use m_general_force_tensor, only : read_Ftensor_input
        
        type(io_H),intent(out)      :: H_io
        character(*),parameter      :: fname='input'
        integer                     :: io_param

        io_param=open_file_read(fname)

        Call read_anisotropy_input(io_param,fname,H_io%aniso)
        Call read_zeeman_input(io_param,fname,H_io%zeeman)
        Call read_J_input(io_param,fname,H_io%J)
        Call read_D_input(io_param,fname,H_io%D)
        Call read_TJ_input(io_param,fname,H_io%TJ)
        Call read_ME_J_input(io_param,fname,H_io%ME_J)
        Call read_ME_D_input(io_param,fname,H_io%ME_D)
        Call read_F_input(io_param,fname,H_io%F)
        Call read_ASR_Ph_input(io_param,fname,H_io%ASR_ph)
        Call read_stark_input(io_param,fname,H_io%stark)
        call read_Mag_Biq_input(io_param,fname,H_io%M_biq)
        call read_sp4_input(io_param,fname,H_io%sp4)
        call read_dip_input(io_param,fname,H_io%dip)
        call read_ExchG_input(io_param,fname,H_io%Exchten)
        call read_PH4_input(io_param,fname,H_io%Ph4)
        call read_SC_input(io_param,fname,H_io%SC)
        call read_Ph_Biq_input(io_param,fname,H_io%U_biq)
        call read_Ftensor_input(io_param,fname,H_io%U_foten)

        call close_file(fname,io_param)
    end subroutine


end module
