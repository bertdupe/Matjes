module m_excitation_shape
use m_excitation_shape_base, only: excitation_shape
use m_type_lattice, only: dim_modes_inner,op_name_to_int
use m_heaviside,only:  heaviside_read_string
use m_rampe    ,only:      rampe_read_string
use m_EMwave   ,only:     EMwave_read_string
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit

private
public  ::  read_excitation_shape, excitation_shape

contains

subroutine read_excitation_shape(io,fname,shapes)
    use m_io_read_util
    integer,intent(in)                                  :: io
    character(len=*),intent(in)                         :: fname
    type(excitation_shape),intent(inout),allocatable    :: shapes(:)

    character(len=*),parameter  :: var_name='excitation_shape_t'
    character(len=1000)         :: str
    logical                     :: success
    integer                     :: nread
    integer                     :: stat
    type(excitation_shape)      :: io_exc
    integer                     :: i


    inquire(io,opened=success)
    if(.not.success)then
        write(error_unit,'(A)') "Error reading the excitation shape"
        write(error_unit,'(A)') "Input io-unit is not opened"
        STOP
    endif

    Call set_pos_entry(io,fname,var_name,success)
    if(.not.success)then
        write(output_unit,'(//A)') "No excitation_shape provided in input"
        write(output_unit,'(A/)') "Need to set excitation_shape, if excitations are used"
        return
    endif
    read(io,'(a)',iostat=stat)! nothing to read in "excitation_norm" line, one could put the number there
    nread=0
    do 
        read(io,'(a)',iostat=stat) str
        if (stat < 0)then
            write(error_unit,'(A)') "end of input file reached reading excitation"
            exit
        endif
        Call shape_read_string(io_exc,str,success)
        if(success)then
            nread=nread+1   
            Call io_exc%unset()
        else
            exit
        endif
    enddo

    if(nread>0)then
        write(output_unit,'(/AI3A/)') "Found ",nread," excitation_shape entries which are read now"
        !return do beginning of excitation data
        do i=1,Nread+1
            backspace(io)
        enddo
        allocate(shapes(nread))
        do i=1,Nread
            read(io,'(a)',iostat=stat) str
            Call shape_read_string(shapes(i),str,success)
            if(.not.success) ERROR STOP "PROGRAMMING MISTAKE, THIS SHOULD ALWAYS WORK"
        enddo
    else
        write(error_unit,'(2/A/A/)') "WARNING, specified excitation_shape, but no excitation_norms are found", "CHECK INPUT"
    endif

end subroutine

subroutine shape_read_string(this,string,success)
    class(excitation_shape),intent(inout)   :: this
    character(len=*),intent(in)             :: string
    logical,intent(out)                     :: success

    integer     ::  stat

    character(len=100)                      :: dummy_name
    character(len=100)                      :: shape_name
    integer                                 :: size_value
    integer                                 :: pos
    integer                                 :: Nreal

    success=.false.
    read(string,*,iostat=stat)  dummy_name, shape_name
    if(stat/=0) return

    this%name=trim(dummy_name)

    select case(trim(shape_name))
    case('heaviside')
        Call heaviside_read_string(this,string)
    case('rampe')
        Call rampe_read_string(this,string)
    case('EMwave')
        Call Emwave_read_string(this,string)
    case default
        write(error_unit,'(A)') "Error reading excitation shape"
        write(error_unit,'(2A)') "Shape identifier not implemented:",trim(shape_name)
        write(error_unit,'(A)') "Error reading line (second entry):"
        write(error_unit,'(A)') string
        STOP
    end select
    if(stat/=0)then
        write(error_unit,'(A)') "Error reading excitation norm"
        write(error_unit,'(A)') "Failed to read all information from line:"
        write(error_unit,'(A)') string
        write(error_unit,'(A,I3,A)') "Two strings should be followed by ",Nreal," reals, which are not recognized"
        STOP
    endif
    success=.true.
end subroutine
end module
