module m_excitation_io
use m_type_lattice, only: order_parameter_name, op_name_to_int
use m_io_read_util
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit

private
public excitation_io, read_excitation_io, check_excitation_io
type excitation_io
    integer             :: op   !operator index
    character(len=100)  :: shape_t_name
    character(len=100)  :: shape_r_name
contains
    procedure   ::   read_string
end type

character(len=*),parameter  :: var_name='excitations'

contains

subroutine check_excitation_io(io,fname,found)
    integer,intent(in)              :: io
    character(len=*),intent(in)     :: fname
    logical,intent(out)             :: found 
    inquire(io,opened=found)
    if(.not.found)then
        write(error_unit,'(A)') "Error reading the excitations"
        write(error_unit,'(A)') "Input io-unit is not opened"
        STOP
    endif

    Call set_pos_entry(io,fname,var_name,found)
    if(.not.found) write(output_unit,'(/A/)') "No excitations provided in input"
end subroutine
subroutine read_excitation_io(io,fname,excitations)
    integer,intent(in)                              :: io
    character(len=*),intent(in)                     :: fname
    type(excitation_io),intent(inout),allocatable   :: excitations(:)

    character(len=1000)         :: str
    logical                     :: success
    integer                     :: nread
    integer                     :: stat
    type(excitation_io)         :: io_exc
    integer                     :: i


    inquire(io,opened=success)
    if(.not.success)then
        write(error_unit,'(A)') "Error reading the excitations"
        write(error_unit,'(A)') "Input io-unit is not opened"
        STOP
    endif

    Call set_pos_entry(io,fname,var_name,success)
    if(.not.success)then
        write(output_unit,'(/A/)') "No excitations provided in input"
        return
    endif
    read(io,'(a)',iostat=stat)! nothing to read in "excitations" line, one could put the number there

    nread=0
    do 
        read(io,'(a)',iostat=stat) str
        if (stat < 0)then
            write(error_unit,'(A)') "end of input file reached reading excitation"
            exit
        endif
        Call read_string(io_exc,trim(str),success)
        if(success)then
            nread=nread+1   
        else
            exit
        endif
    enddo

    if(nread>0)then
        write(output_unit,'(/A,I3,A/)') "Found ",nread," excitation entries which are read now"
        !return do beginning of excitation data
        do i=1,Nread+1
            backspace(io)
        enddo
        allocate(excitations(nread))
        do i=1,Nread
            read(io,'(a)',iostat=stat) str
            Call read_string(excitations(i),trim(str),success)
            if(.not.success) ERROR STOP "PROGRAMMING MISTAKE, THIS SHOULD ALWAYS WORK"
        enddo
    else
        write(error_unit,'(2/A/A/)') "WARNING, specified excitations, but no excitations are found", "CHECK INPUT"
    endif
end subroutine

subroutine read_string(this,string,success)
    class(excitation_io),intent(inout)      :: this
    character(len=*),intent(in)             :: string
    logical,intent(out)                     :: success

    integer             ::  stat
    character(len=100)  :: operator_name

    success=.false.
    read(string,*,iostat=stat)  operator_name, this%shape_t_name, this%shape_r_name
    if(stat/=0) return
    
    if(.not.any(operator_name==order_parameter_name))then
        write(error_unit,'(/2A)') "ERROR, Found excitation operator name:",trim(operator_name)
        write(error_unit,'(A)') "This operator is not defined in the lattice type"
        write(error_unit,'(A)') "Error appeared in line:"
        write(error_unit,'(A)') string
    else
        this%op=op_name_to_int(trim(operator_name))
        success=.true.
    endif
end subroutine
end module
