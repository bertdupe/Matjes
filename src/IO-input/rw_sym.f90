module m_rw_sym
implicit none
contains
subroutine rw_sym(mode,cal_sym,tol_sym)
    use m_derived_types
    use m_io_files_utils
    use m_io_utils
    integer, intent(inout) :: mode
    logical, intent(inout) :: cal_sym
    real(8), intent(inout) :: tol_sym
    ! internal
    integer :: io_input

    mode=1
    cal_sym=.true.
    io_input=open_file_read('input')
    call get_parameter(io_input,'input','cal_sym',cal_sym)
    call get_parameter(io_input,'input','sym_mode',mode)
    call get_parameter(io_input,'input','tol_sym',tol_sym)
    call close_file('input',io_input)

end subroutine rw_sym


end module
