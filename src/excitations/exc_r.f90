module m_exc_r
use m_type_lattice, only: dim_modes_inner,op_name_to_int
use m_constants
use m_laguerre_interface
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit

private
public  ::  read_excitation_shape_r, excitation_shape_r

type excitation_shape_r
    character(len=30)                    :: name='plane' !name for identification     
    real(8)                              :: center(3)=0.0d0
    real(8)                              :: cutoff=1.0d0
    integer                              :: l=0, p=0 ! for Laguerre-Gauss
    ! integer                              :: m=0, n=0 ! for Hermite-Gauss
    real(8)                              :: width_0, wavelength ! for Laguerre-Gauss, width of the beam at its waist and wavelength
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
        real(8)                                 :: shape_r(2) ! complex number (magnitude and phase)
                                                              ! shape_r(1) * exp(i*shape_r(2))
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
        write(output_unit,'(/A,I3,A/)') "Found ",nread," excitation_shape_r entries which are read now"
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
    real(8)                                 :: period

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
    case('laguerre-gauss')
        Nreal=7
        read(string,*,iostat=stat)  dummy_name, shape_r_name, this%center, this%l, this%p, this%width_0, period
        this%wavelength = c * period
        this%shape_r => shape_r_laguerre_gauss
        this%print_r => print_r_laguerre_gauss
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
    real(8)                                 :: shape_r(2)
    shape_r(1)=1.0d0
    shape_r(2)=0.0d0
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
    real(8)                                 :: shape_r(2)

    shape_r(1)=0.0d0
    if (all(abs(R-this%center).lt.this%cutoff)) shape_r(1)=1.0d0
    shape_r(2)=0.0d0
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
    real(8)                                 :: shape_r(2)
!    real(8)     :: dist
    
!    dist=shape_r2(R-R0)
!    shape_r=0.5d0*(sign(1.0d0,cutoff-dist)+1.0d0)
    shape_r(1)=0.0d0
    if (norm2(R-this%center).lt.this%cutoff) shape_r(1)=1.0d0
    shape_r(2)=0.0d0
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
    real(8)                                 :: shape_r(2)
    real(8)     :: tmp

    tmp=norm2(R-this%center)
    tmp=tmp/this%cutoff
    tmp=-tmp*tmp
    tmp=max(tmp,-200.0d0)  !prevent exp(tmp) underflow 
    shape_r(1)=exp(tmp)
    shape_r(2)=0.0d0
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

function shape_r_laguerre_gauss(this,R)result(shape_r)
    class(excitation_shape_r),intent(in)    :: this
    real(8), intent(in)                     :: R(3)
    real(8)                                 :: shape_r(2)
    real(8)     :: tmp
    real(8)     :: radius, phi, z ! cylindrical coordinates
    real(8)     :: w_z    ! w_z: width of beam at z
    real(8)     :: R_curv ! radius of curvature
    real(8)     :: z_R    ! Rayleigh length
    real(8)     :: psi    ! Gouy phase
    real(8)     :: normalization_factor
    real(8)     :: laguerre
    real(8)     :: magnitude, phase ! magn. and phase of shape_r

    z = R(3) - this%center(3)
    radius = norm2(R(1:2)-this%center(1:2))
    phi = datan2(R(2)-this%center(2), R(1)-this%center(1))
    z_R = pi * this%width_0**2 / this%wavelength
    w_z = this%width_0 * sqrt(1.0d0 + (z/z_R)**2)
    if(z /= 0.0d0) R_curv = z * (1.0d0 + (z_R/z)**2)
    psi = atan(z/z_R)
    normalization_factor = sqrt(2**(abs(this%l)+1.0d0) * gamma(real(this%p+1.0d0)) / (pi * gamma(real(this%p+abs(this%l)+1.0d0))))
    Call laguerre_polynomial_scalar(abs(this%l), this%p, 2.0d0*radius**2/w_z**2, laguerre)
    ! print *, "2.0d0*radius**2/w_z**2: ", 2.0d0*radius**2/w_z**2
    ! print *, "laguerre: ", laguerre
    if(radius**2/w_z**2 >= 200.0d0) then ! to prevent from exponential underflow
        magnitude = 0
        phase = 0
    else
        magnitude = normalization_factor/w_z * (radius*sqrt(2.0d0)/w_z)**abs(this%l) * laguerre * exp(-radius**2/w_z**2)
        if(z /= 0.0d0) then
            phase = pi*radius**2/(this%wavelength*R_curv) + 2*pi*z/this%wavelength + this%l*phi - (2*this%p+abs(this%l)+1)*psi
        else
            phase = 2*pi*z/this%wavelength + this%l*phi - (2*this%p+abs(this%l)+1)*psi
        endif
    endif
    shape_r(1) = magnitude
    shape_r(2) = phase
end function

subroutine print_r_laguerre_gauss(this,io)
    class(excitation_shape_r),intent(in)    :: this
    integer,intent(in)                      :: io

    write(io,'(3X,A)') "Real-space shape: laguerre-gauss (propagating in the +z direction)"
    write(io,'(6X,A)') "Equation: ..."
    write(io,'(6X,A)') "Parameters:"
    write(io,'(9X,A,I3,A,I3)')  "mode                  : l = ",this%l,", p = ", this%p
    write(io,'(9X,A,3F16.8,A)') "center position       : ",this%center," nm"
    write(io,'(9X,A,F16.8,A)')  "width (w_0)           : ",this%width_0," nm"
    write(io,'(9X,A,F16.8,A)')  "Wavelength (lambda) : ",this%wavelength," nm"
end subroutine

end module
