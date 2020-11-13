module m_rw_H
    use m_input_H_types
    public

    contains
    subroutine rw_H(H_io)
        use m_io_files_utils, only : open_file_read,close_file
        use m_anisotropy_heisenberg, only: read_anisotropy_input
        use m_zeeman, only: read_zeeman_input
        use m_couplage_ME, only: read_ME_input
        use m_exchange_heisenberg_J, only: read_J_input
        use m_exchange_heisenberg_D, only: read_D_input
        
        type(io_H),intent(out)      :: H_io
        character(*),parameter      :: fname='input'
        integer                     :: io_param

        io_param=open_file_read(fname)

        Call read_anisotropy_input(io_param,fname,H_io%aniso)
        Call read_zeeman_input(io_param,fname,H_io%zeeman)
        Call read_ME_input(io_param,fname,H_io%ME)
        Call read_J_input(io_param,fname,H_io%J)
        Call read_D_input(io_param,fname,H_io%D)

        call close_file(fname,io_param)


    end subroutine


end module
