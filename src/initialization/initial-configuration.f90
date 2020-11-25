module m_init_config
use m_derived_types, only: lattice
use m_io_utils,only: get_parameter
use m_input_types, only: extpar_input
private
public :: init_config_lattice
contains

subroutine init_config_lattice(lat,initialized,extpar_io,fname_in)
    !loops through all order parameters and checks if there is some initial configuration supplied from the file input
    use m_io_files_utils,only : open_file_read,close_file
    use m_type_lattice,only: number_different_order_parameters,order_parameter_name
    type(lattice), intent(inout)        :: lat
    logical,intent(inout)               :: initialized(number_different_order_parameters)
    type(extpar_input),intent(in)       :: extpar_io
    character(*), intent(in),optional   :: fname_in
    integer :: io

    character(*),parameter              :: fname_default='init.config'
    character(:), allocatable           :: fname
    logical :: exists
    integer :: i
    integer :: dim_mode
    real(8),pointer ::  state(:)

    if(present(fname_in))then
        fname=fname_in
    else
        fname=fname_default
    endif
    inquire(file=fname,exist=exists)
    if(.not.exists)then
        write(*,'(///A)') 'FILE:',fname," NOT FOUND"
        write(*,'(A///)') 'SKIPPING init_config'
        return
    endif
    io=open_file_read(fname)
    do i=1,number_different_order_parameters
        if(initialized(i)) cycle    !don't try to intialize, if already initialized before
        dim_mode=lat%get_order_dim(i,ignore=.true.)
        if(dim_mode<1) cycle
        Call lat%set_order_point(i,state)
        Call init_config_order(io,fname,lat,trim(order_parameter_name(i)),dim_mode,state,extpar_io,initialized(i))
        nullify(state)
    enddo
    call close_file(fname,io)
end subroutine

subroutine init_config_order(io,fname,lat,ordname,dim_mode,state,extpar_io,init)
    use m_init_default
    use m_init_DW
    use m_init_heavyside
    use m_init_Sk
    use m_init_sky_lin
    use m_init_Sklattice
    use m_init_spiral
    use m_init_random
    integer,intent(in)              :: io       !init-file io-unit
    character(*),intent(in)         :: fname    !init-file name 
    type(lattice), intent(in)       :: lat      !entire lattice containing geometric information
    character(*),intent(in)         :: ordname  !name of the order parameter
    integer,intent(in)              :: dim_mode !dimension of the order parameter in each cell
    real(8),pointer,intent(inout)   :: state(:) !pointer the the order parameter
    type(extpar_input),intent(in)   :: extpar_io
    logical,intent(out)             :: init     !returns if initialization has occured

    character(30)                   :: configuration

    configuration=""
    call get_parameter(io,fname,'configuration_'//ordname,configuration)

    init=.true.
    select case (adjustl(configuration))
        case('spiral')
            call init_spiral(io,fname,lat,ordname,dim_mode,state)
        case('domainwall')
            call init_DW(io,fname,lat,ordname,dim_mode,state)
        case('heavyside')
            call init_heavyside(lat,dim_mode,state)
        case('skyrmion')
            call init_spiral(io,fname,lat,ordname,dim_mode,state)
            call init_Sk(io,fname,lat,ordname,dim_mode,state)
        case('skyrmion_lin')
            Call init_Sky_lin(io,fname,lat,ordname,dim_mode,state)
        case('skyrmionla')
            call init_spiral(io,fname,lat,ordname,dim_mode,state)
            call init_Sk_lattice(io,fname,lat,ordname,dim_mode,state)
        case('random')
            write(6,'(3a)') 'random configuration for ',ordname,' was chosen'
            call init_random(dim_mode,state)
        case('')
            write(6,'(3a)') 'Using default initial configuration for: ',ordname
            call init_default(ordname,dim_mode,state,extpar_io)
        case default
            write(6,'(3a)') 'Unknown initial configuration for ',ordname,' was found'
            write(6,'(2a)') 'initial configuration: ', adjustl(configuration)
            init=.false.
    end select
end subroutine 

end module m_init_config
