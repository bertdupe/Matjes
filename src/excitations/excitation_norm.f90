module m_excitation_norm
use m_type_lattice, only: dim_modes_inner,op_name_to_int
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit

private
public  ::  read_excitation_norms, excitation_norm

type excitation_norm
    character(len=30)                    :: name='plane' !name for identification     
    real(8)                              :: center(3)=0.0d0
    real(8)                              :: cutoff=1.0d0
    procedure(int_norm), pointer,pass    :: norm=>norm_plane
contains
    procedure   ::   read_string
end type

abstract interface
    function int_norm(this,R)result(norm)
        import excitation_norm
        class(excitation_norm),intent(in)   :: this
        real(8), intent(in)                 :: R(3)
        real(8)                             :: norm
    end function
end interface

contains

subroutine read_excitation_norms(io,fname,norm)
    use m_io_read_util
    integer,intent(in)                              :: io
    character(len=*),intent(in)                     :: fname
    type(excitation_norm),intent(inout),allocatable :: norm(:)

    character(len=*),parameter  :: var_name='excitation_shape_r'
    character(len=100)          :: str
    logical                     :: success
    integer                     :: nread
    integer                     :: stat
    type(excitation_norm)       :: io_exc
    integer                     :: i


    inquire(io,opened=success)
    if(.not.success)then
        write(error_unit,'(A)') "Error reading the excitation norms"
        write(error_unit,'(A)') "Input io-unit is not opened"
        STOP
    endif

    Call set_pos_entry(io,fname,var_name,success)
    if(.not.success)then
        write(output_unit,'(/A/)') "No excitation_norm provided in input"
        write(output_unit,'(A)') "Setting only plane as default"
        allocate(norm(1))
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
        Call io_exc%read_string(str,success)
        if(success)then
            nread=nread+1   
        else
            exit
        endif
    enddo

    if(nread>0)then
        write(output_unit,'(/AI3A/)') "Found ",nread," excitation_norm entries which are read now"
        !return do beginning of excitation data
        do i=1,Nread+1
            backspace(io)
        enddo
        allocate(norm(nread))
        do i=1,Nread
            read(io,'(a)',iostat=stat) str
            Call norm(i)%read_string(str,success)
            if(.not.success) ERROR STOP "PROGRAMMING MISTAKE, THIS SHOULD ALWAYS WORK"
        enddo
    else
        write(error_unit,'(2/A/A/)') "WARNING, specified excitation_norm, but no excitation_norms are found", "CHECK INPUT"
    endif

end subroutine

subroutine read_string(this,string,success)
    class(excitation_norm),intent(inout)    :: this
    character(len=*),intent(in)             :: string
    logical,intent(out)                     :: success

    integer     ::  stat

    character(len=100)                      :: dummy_name
    character(len=100)                      :: norm_name
    integer                                 :: size_value
    integer                                 :: pos
    integer                                 :: Nreal

    success=.false.
    read(string,*,iostat=stat)  dummy_name, norm_name
    if(stat/=0) return

    this%name=trim(dummy_name)

    select case(trim(norm_name))
    case('plane')
        !no additional data has to be read for plane
        Nreal=0
        this%norm => norm_plane
    case('square')
        Nreal=4
        read(string,*,iostat=stat)  dummy_name, norm_name, this%center,this%cutoff
        this%norm => norm_square
    case('cylinder')
        Nreal=4
        read(string,*,iostat=stat)  dummy_name, norm_name, this%center,this%cutoff
        this%norm => norm_cylinder
    case('gaussian')
        Nreal=4
        read(string,*,iostat=stat)  dummy_name, norm_name, this%center,this%cutoff
        this%norm => norm_gaussian
    case default
        write(error_unit,'(A)') "Error reading excitation norm"
        write(error_unit,'(2A)') "Norm identifier not implemented (second entry):", trim(norm_name)
        write(error_unit,'(A)') "Error reading line:"
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


function norm_plane(this,R)result(norm)
    class(excitation_norm),intent(in)   :: this
    real(8), intent(in)                 :: R(3)
    real(8)                             :: norm
    norm=1.0d0

end function

function norm_square(this,R)result(norm)
    class(excitation_norm),intent(in)   ::  this
    real(8), intent(in)                 :: R(3)
    real(8)                             :: norm

    norm=0.0d0
    if (all(abs(R-this%center).lt.this%cutoff)) norm=1.0d0
end function

function norm_cylinder(this,R)result(norm)
    !actually a sphere, but I will not change the functionality at this point
    class(excitation_norm),intent(in)   :: this
    real(8), intent(in)                 :: R(3)
    real(8)                             :: norm
!    real(8)     :: dist
    
!    dist=norm2(R-R0)
!    norm=0.5d0*(sign(1.0d0,cutoff-dist)+1.0d0)
    norm=0.0d0
    if (norm2(R-this%center).lt.this%cutoff) norm=1.0d0
end function

function norm_gaussian(this,R)result(norm)
    class(excitation_norm),intent(in)   :: this
    real(8), intent(in)                 :: R(3)
    real(8)                             :: norm
    real(8)     :: tmp

    tmp=norm2(R-this%center)
    tmp=tmp/this%cutoff
    tmp=-tmp*tmp
    tmp=max(tmp,-200.0d0)  !prevent exp(tmp) underflow 
    norm=exp(tmp)
end function
end module
