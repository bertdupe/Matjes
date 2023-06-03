module m_io_read_util
!module with some general routines to reduce code duplication in file reading
!so far mostly used in TB multiple readin
public

contains

function is_comment(line)
    !Function which checks if a line starts with a comment
    character(*),intent(in) :: line             !Line, that is checked if it is a comment
    logical                 :: is_comment

    character(len=100)  :: test_comment, dummy  !temporary chars for reading
    integer             :: stat

    !Array that contains all chars which indicate a line as a comment
    character(len=1),parameter  ::  comment_chars(1) = ["#"]    

    read(line,*,iostat=stat) test_comment, dummy
    is_comment=any(test_comment(1:1) .eq. comment_chars)
end function

subroutine set_pos_entry(io,fname,var_name,success_out)
    !sets the io-unit (io) to the line starting with var_name
    !If success_out is supplied, it return whether such a line is found
    !and sets the io-unit to that line or the end
    !If success_out is not present this crashes if this line is not found
    use, intrinsic :: iso_fortran_env, only : error_unit
    integer, intent(in)                         :: io
    character(len=*), intent(in)                :: fname,var_name
    logical,intent(out),optional                :: success_out
    logical     ::  success
    character(len=100) :: str
    integer :: length_string, stat

    length_string=len_trim(var_name)
    success=.false.
    Call check_io_fname(io,fname,var_name)

    rewind(io)
    do
        read (io,'(a)',iostat=stat) str
        if (stat /= 0) exit
        str= trim(adjustl(str))
        if (len_trim(str)==0) cycle
        if (str(1:1) == '#' ) cycle
        if (len_trim(str)<length_string) cycle
        if (len_trim(str)>length_string)then
            if(str(length_string+1:length_string+1)/=" ") cycle
        endif

        if (str(1:length_string) == var_name(1:length_string))then
            success=.true.
            backspace(io)
            exit
        endif
    enddo
    if(present(success_out))then
        success_out=success
    else if(.not.success)then
        write(error_unit,'(/2A)') "failed to find necessary keyword: ", var_name
        STOP "ERROR, Fix input?"
    endif
end subroutine

subroutine check_io_fname(io,fname,var_name)
    !checks if the io-unit (io) is indeed associated with fname, if not this will throw an error
    use, intrinsic :: iso_fortran_env, only : error_unit
    integer, intent(in)             :: io
    character(len=*), intent(in)    :: fname,var_name

    integer :: test_unit
    logical :: isopen
    inquire(file=fname, number=test_unit, opened=isopen)
    if(.not.isopen.or.test_unit/=io)then
        write(error_unit,*) 'io-unit and fname: "',fname,'" do not match'
        write(error_unit,*) 'error occured reading var_name: "',var_name,'"'
        ERROR STOP "Error in input or programming mistake"
    endif
end subroutine

subroutine check_further_entry(io,fname,var_name,another_entry_out)
    !checks if there is another entry of var_name in the io io-unit
    !if another_entry_out is supplied, the result is return, otherwise
    !this throws an error if another entry is found
    use, intrinsic :: iso_fortran_env, only : error_unit
    integer, intent(in)                         :: io
    character(len=*), intent(in)                :: fname,var_name
    logical,intent(out),optional                :: another_entry_out  !true if another entry is found, omit for crash

    logical            :: another_entry
    character(len=100) :: str
    integer :: length_string, stat,test_unit
    logical :: isopen

    length_string=len_trim(var_name)
    another_entry=.false.
    inquire(file=fname, number=test_unit, opened=isopen)
    if(.not.isopen.or.test_unit/=io) ERROR STOP "IO-unit is not open or wrongly assigned"
    do
        read (io,'(a)',iostat=stat) str
        if (stat /= 0) exit
        str= trim(adjustl(str))
        if (len_trim(str)==0) cycle
        if (str(1:1) == '#' ) cycle
        if (len_trim(str)<length_string) cycle
        if (len_trim(str)>length_string)then
            if(str(length_string+1:length_string+1)/=" ") cycle
        endif

        if (str(1:length_string) == var_name(1:length_string))then
            another_entry=.true.
            exit
        endif
    enddo
    if(another_entry)then
        if(present(another_entry_out))then
            another_entry_out=another_entry
        else
            write(error_unit,'(/2A)') "ERROR, found unexpected second entry of: ",var_name
            STOP "ABORTING, INPUT INCORRECT?"
        endif
    endif
end subroutine

subroutine write_info_number_found(Nentry,Nnonzero,var_name)
    use, intrinsic :: iso_fortran_env, only : error_unit,output_unit
    !write info how many entries have been found
    integer,intent(in)  ::  Nentry,Nnonzero
    character(len=*)    ::  var_name
    if(Nentry<1)then
        write(error_unit,'(/2A/A/)') "Found no entries for ",var_name,' although the keyword is specified'
        ERROR STOP "INPUT PROBABLY WRONG"
    endif
    if(Nnonzero<1)then
        write(error_unit,'(/2A/A/)') "Found no nonzero entries for ",var_name,' although the keyword is specified'
        ERROR STOP "INPUT PROBABLY WRONG"
    endif
    write(output_unit,'(/A,I6,2A)') "Found ",Nnonzero," nonzero entries for input-parameter: ",var_name
end subroutine

end module
