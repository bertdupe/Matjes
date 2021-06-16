module m_exc_r
use m_type_lattice, only: dim_modes_inner,op_name_to_int
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit

private
public  ::  read_excitation_shape_r, excitation_shape_r

type excitation_shape_r
    character(len=30)                    :: name='plane' !name for identification     
    real(8)                              :: center(3)=0.0d0
    real(8)                              :: cutoff=1.0d0
    procedure(int_shape_r), pointer,pass :: shape_r=>shape_r_plane
    procedure(int_print_r), pointer,pass :: print_r=>print_r_plane
contains
    procedure   ::   read_string
end type

abstract interface
    function int_shape_r(this,R)result(shape_r)
        import excitation_shape_r
        class(excitation_shape_r),intent(in)    :: this
        real(8), intent(in)                     :: R(3)
        real(8)                                 :: shape_r
    end function
    subroutine int_print_r(this,io)
        import excitation_shape_r
        class(excitation_shape_r),intent(in)    :: this
        integer,intent(in)                      :: io
    end subroutine 
end interface

contains

subroutine read_excitation_shape_r(io,fname,shape_r)
    use m_io_read_util
    integer,intent(in)                              :: io
    character(len=*),intent(in)                     :: fname
    type(excitation_shape_r),intent(inout),allocatable :: shape_r(:)

    character(len=*),parameter  :: var_name='excitation_shape_r'
    character(len=100)          :: str
    logical                     :: success
    integer                     :: nread
    integer                     :: stat
    type(excitation_shape_r)       :: io_exc
    integer                     :: i


    inquire(io,opened=success)
    if(.not.success)then
        write(error_unit,'(A)') "Error reading the excitation shape_r"
        write(error_unit,'(A)') "Input io-unit is not opened"
        STOP
    endif

    Call set_pos_entry(io,fname,var_name,success)
    if(.not.success)then
        write(output_unit,'(/A/)') "No excitation_shape_r provided in input"
        write(output_unit,'(A)') "Setting only plane as default"
        allocate(shape_r(1))
        return
    endif
    read(io,'(a)',iostat=stat)! nothing to read in "excitation_shape_r" line, one could put the number there

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
        write(output_unit,'(/AI3A/)') "Found ",nread," excitation_shape_r entries which are read now"
        !return do beginning of excitation data
        do i=1,Nread+1
            backspace(io)
        enddo
        allocate(shape_r(nread))
        do i=1,Nread
            read(io,'(a)',iostat=stat) str
            Call shape_r(i)%read_string(str,success)
            if(.not.success) ERROR STOP "PROGRAMMING MISTAKE, THIS SHOULD ALWAYS WORK"
        enddo
    else
        write(error_unit,'(2/A/A/)') "WARNING, specified excitation_shape_r, but no excitation_shape_r are found", "CHECK INPUT"
    endif

end subroutine

subroutine read_string(this,string,success)
    class(excitation_shape_r),intent(inout) :: this
    character(len=*),intent(in)             :: string
    logical,intent(out)                     :: success

    integer     ::  stat

    character(len=100)                      :: dummy_name
    character(len=100)                      :: shape_r_name
    integer                                 :: size_value
    integer                                 :: pos
    integer                                 :: Nreal

    success=.false.
    read(string,*,iostat=stat)  dummy_name, shape_r_name
    if(stat/=0) return

    this%name=trim(dummy_name)

    select case(trim(shape_r_name))
    case('plane')
        !no additional data has to be read for plane
        Nreal=0
        this%shape_r => shape_r_plane
        this%print_r => print_r_plane
    case('square')
        Nreal=4
        read(string,*,iostat=stat)  dummy_name, shape_r_name, this%center,this%cutoff
        this%shape_r => shape_r_square
        this%print_r => print_r_square
    case('cylinder')
        Nreal=4
        read(string,*,iostat=stat)  dummy_name, shape_r_name, this%center,this%cutoff
        this%shape_r => shape_r_sphere
        this%print_r => print_r_sphere
    case('gaussian')
        Nreal=4
        read(string,*,iostat=stat)  dummy_name, shape_r_name, this%center,this%cutoff
        this%shape_r => shape_r_gaussian
        this%print_r => print_r_gaussian
    case default
        write(error_unit,'(A)') "Error reading excitation_shape_r"
        write(error_unit,'(2A)') "shape_r identifier not implemented (second entry):", trim(shape_r_name)
        write(error_unit,'(A)') "Error reading line:"
        write(error_unit,'(A)') string
        STOP
    end select
    if(stat/=0)then
        write(error_unit,'(A)') "Error reading excitation_shape_r"
        write(error_unit,'(A)') "Failed to read all information from line:"
        write(error_unit,'(A)') string
        write(error_unit,'(A,I3,A)') "Two strings should be followed by ",Nreal," reals, which are not recognized"
        STOP
    endif
    success=.true.
end subroutine


function shape_r_plane(this,R)result(shape_r)
    class(excitation_shape_r),intent(in)    :: this
    real(8), intent(in)                     :: R(3)
    real(8)                                 :: shape_r
    shape_r=1.0d0
end function

subroutine print_r_plane(this,io)
    class(excitation_shape_r),intent(in)    :: this
    integer,intent(in)                      :: io

    write(io,'(3X,A)') "Real-space shape: plane"
    write(io,'(6X,A)') "No necessary parameters"
end subroutine
    

function shape_r_square(this,R)result(shape_r)
    class(excitation_shape_r),intent(in)    :: this
    real(8), intent(in)                     :: R(3)
    real(8)                                 :: shape_r

    shape_r=0.0d0
    if (all(abs(R-this%center).lt.this%cutoff)) shape_r=1.0d0
end function

subroutine print_r_square(this,io)
    class(excitation_shape_r),intent(in)    :: this
    integer,intent(in)                      :: io

    write(io,'(3X,A)') "Real-space shape: square"
    write(io,'(6X,A)') "Parameters:"
    write(io,'(9X,A,3F16.8,A)') "center position: ",this%center," nm"
    write(io,'(9X,A,F16.8,A)')  "half width     : ",this%cutoff," nm"
end subroutine


function shape_r_sphere(this,R)result(shape_r)
    !actually a sphere, but I will not change the functionality at this point
    class(excitation_shape_r),intent(in)    :: this
    real(8), intent(in)                     :: R(3)
    real(8)                                 :: shape_r
!    real(8)     :: dist
    
!    dist=shape_r2(R-R0)
!    shape_r=0.5d0*(sign(1.0d0,cutoff-dist)+1.0d0)
    shape_r=0.0d0
    if (norm2(R-this%center).lt.this%cutoff) shape_r=1.0d0
end function

subroutine print_r_sphere(this,io)
    class(excitation_shape_r),intent(in)    :: this
    integer,intent(in)                      :: io

    write(io,'(3X,A)') "Real-space shape: sphere (called cylinder)"
    write(io,'(6X,A)') "Parameters:"
    write(io,'(9X,A,3F16.8,A)') "center position: ",this%center," nm"
    write(io,'(9X,A,F16.8,A)')  "radius         : ",this%cutoff," nm"
end subroutine


function shape_r_gaussian(this,R)result(shape_r)
    class(excitation_shape_r),intent(in)    :: this
    real(8), intent(in)                     :: R(3)
    real(8)                                 :: shape_r
    real(8)     :: tmp

    tmp=norm2(R-this%center)
    tmp=tmp/this%cutoff
    tmp=-tmp*tmp
    tmp=max(tmp,-200.0d0)  !prevent exp(tmp) underflow 
    shape_r=exp(tmp)
end function

subroutine print_r_gaussian(this,io)
    class(excitation_shape_r),intent(in)    :: this
    integer,intent(in)                      :: io

    write(io,'(3X,A)') "Real-space shape: gaussian"
    write(io,'(6X,A)') "Equation: e^(-((pos-center)/w)^2)"
    write(io,'(6X,A)') "Parameters:"
    write(io,'(9X,A,3F16.8,A)') "center position: ",this%center," nm"
    write(io,'(9X,A,F16.8,A)')  "width (w)      : ",this%cutoff," nm"
end subroutine

end module
