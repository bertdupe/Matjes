module m_excitations
use m_simu_parameters, only : type_excitations
!use m_shape_excitations
use m_get_position
use m_convert
use m_type_lattice
use m_exc_r
use m_exc_t
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit
implicit none

! variable that contains the the excitations form (sweep of EM field...)
type excitations
    real(8)           :: start_value(maxval(dim_modes_inner))
    real(8)           :: end_value(maxval(dim_modes_inner))
    integer           :: size_value
    real(kind=8)      :: t_start,t_end
! name of the order parameter to change
    character(len=30) :: name
    integer           :: op !order parameter as ordered by m_type_lattice
end type excitations

type excitation
    integer     :: op=0
    integer     :: int_shape_t=0
    integer     :: int_shape_r=0
end type

type excitation_combined
    type(excitation),allocatable            :: exc(:)
    type(excitation_shape_r),allocatable    :: shape_r(:)
    type(excitation_t),allocatable          :: shape_t(:)

    logical :: op_used(number_different_order_parameters)=.false. 
    real(8), allocatable :: position(:,:)
end type

private
public :: read_all_excitations, update_exc, excitation_combined

contains

subroutine read_all_excitations(fname,excitations)
    use m_io_utils
    use m_io_files_utils
    use m_io_read_util
    use m_excitation_io
    character(len=*),intent(in)             :: fname
    type(excitation_combined),intent(inout) :: excitations

    type(excitation_io),allocatable     :: io_exc(:)

    integer             :: i,j, Nexc, io
   
    open(newunit=io,file=fname)
    !get different excitation shape_r
    Call read_excitation_shape_r(io,fname,excitations%shape_r)
    
    !get different excitation shape_t
    Call read_excitation_shape_t(io,fname,excitations%shape_t)

    !read excitation combinations
    Call read_excitation_io(io,fname,io_exc)

    Nexc=0
    if(allocated(io_exc)) Nexc=size(io_exc) 
    if(Nexc==0)then
        !no excitations inserted
        !clean up and leave
        if(allocated(excitations%shape_t)) deallocate(excitations%shape_t)
        if(allocated(excitations%shape_r)) deallocate(excitations%shape_r)
        return  
    else
        !check allocations
        if(.not.allocated(excitations%shape_t))then
            write(error_unit,'(A)') "ERROR, Found some excitations, but 'excitation_shape_t' is not set"
            write(error_unit,'(A)') "CHECK INPUT"
            STOP
        endif
        if(.not.allocated(excitations%shape_r))then
            write(error_unit,'(A)') "ERROR, Found some excitations, but 'excitation_shape_r' is not set"
            write(error_unit,'(A)') "CHECK INPUT"
            STOP
        endif
    endif

    !set excitation entries
    allocate(excitations%exc(Nexc))
    do i=1,Nexc
        excitations%exc(i)%op=io_exc(i)%op
        !find associated shape entry
        do j=1,size(excitations%shape_t)
            if(io_exc(i)%shape_t_name==excitations%shape_t(j)%name)then
                if(excitations%shape_t(j)%dim_mode/=dim_modes_inner(io_exc(i)%op))then
                    write(error_unit,'(A)') "ERROR, Excitation operator dimension does not fit associated excitation shape"
                    write(error_unit,'(2(AI3))') "Error associating excitation ",i,"with shape number",j
                    STOP
                endif
                excitations%exc(i)%int_shape_t=j
                exit
            endif
        enddo
        if(excitations%exc(i)%int_shape_t==0)then
            write(error_unit,'(A)') "ERROR, Excitation operator shape_t name not found in exitation shape_t"
            write(error_unit,'(2A)') "Failed to find name: ", io_exc(i)%shape_t_name
            write(error_unit,'(A)') "CHECK INPUT"
            STOP
        endif


        !find associated shape_r entry
        do j=1,size(excitations%shape_r)
            if(io_exc(i)%shape_r_name==excitations%shape_r(j)%name)then
                excitations%exc(i)%int_shape_r=j
                exit
            endif
        enddo
        if(excitations%exc(i)%int_shape_r==0)then
            write(error_unit,'(A)') "ERROR, Excitation operator shape_r name not found in exitation shape_r"
            write(error_unit,'(2A)') "Failed to find name: ", io_exc(i)%shape_r_name
            write(error_unit,'(A)') "CHECK INPUT"
            STOP
        endif
    enddo

    do i=1,size(excitations%exc)
        excitations%op_used(excitations%exc%op)=.true.
    enddo

    !terrible position-get
    i=get_lines('positions.dat')
    allocate(excitations%position(3,i),source=0.0d0)
    call get_position(excitations%position,'positions.dat')
end subroutine 

subroutine update_exc(time,lat,dat)
    !routine which updated all order-parameters for which an exciatation in dat is specified
    real(8), intent(in)                     :: time
    type(lattice), intent(inout)            :: lat
    type(excitation_combined),intent(in)    :: dat
    ! internal
    real(8),pointer,contiguous :: opvec(:,:)        !pointer to locally considered mode
    real(8) :: r(3)                                 !local position
    real(8) :: shape_t(maxval(dim_modes_inner))     !prefactor defined by the time (constant in space)
    !help indices
    integer :: op       !index for considered operator as in lattice
    integer :: dim_mode !dimension of considered mode
    integer :: i_r, i_t !integer identifier for the considered shape in real- of time-space 
    integer :: i,j      !loop variables

    !set initial values to zero for considered operators
    do j=1,number_different_order_parameters
        if(dat%op_used(j))then
            Call lat%set_order_point_inner(j,opvec)
            opvec=0.0d0
        endif
    enddo
 
    do j=1,size(dat%exc)
        !set help parameters
        i_r=dat%exc(j)%int_shape_r
        i_t=dat%exc(j)%int_shape_t
        op=dat%exc(j)%op
        dim_mode=dim_modes_inner(op)
        
        !get time-shape 
        shape_t(1:dim_mode)=dat%shape_t(i_t)%get_shape(time)
        if(all(shape_t(1:dim_mode)==0.0d0)) cycle   !nothing to do for this excitation, it is necessarily zero (can to this as soon as it is additive)

        !apply real-space shape
        Call lat%set_order_point_inner(op,opvec)
        do i=1,size(dat%position,2)
            r=dat%position(:,i)
            opvec(:,i)=opvec(:,i)+ shape_t(1:dim_mode) * dat%shape_r(i_r)%shape_r(r)
        enddo
    enddo
    nullify(opvec)
end subroutine
end module m_excitations
