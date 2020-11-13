module m_rw_H
    use m_input_H_types
    public

    contains
    subroutine rw_H(H_io)
        use m_io_files_utils, only : open_file_read,close_file
        use m_anisotropy_heisenberg, only: read_anisotropy_input
        
        type(io_H),intent(out)      :: H_io
        character(*),parameter      :: fname='input'
        integer                     :: io_param

        io_param=open_file_read(fname)

        Call read_anisotropy_input(io_param,fname,H_io%aniso)

        call close_file(fname,io_param)


    end subroutine


end module
